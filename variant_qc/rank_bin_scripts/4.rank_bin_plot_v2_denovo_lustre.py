# Pavlos Antoniou
# 22/01/2021
#  rank and bin variants from RF

from hail import Table
import os
import pprint
from pprint import pformat
import argparse
import hail as hl
import pandas as pd
import numpy as np
import pyspark
from pprint import pformat
import json
import sys
import re
from collections import defaultdict

from pathlib import Path
import logging
from typing import Any, Counter, List, Optional, Tuple, Union, Dict
import uuid
import json
from bokeh.plotting import output_file, save, show
from gnomad.resources.grch38 import gnomad
from gnomad.utils.annotations import unphase_call_expr, add_variant_type
from gnomad.variant_qc.pipeline import create_binned_ht, score_bin_agg
from gnomad.variant_qc.pipeline import test_model, sample_training_examples, get_features_importance

from gnomad.variant_qc.pipeline import train_rf_model
# from gnomad.variant_qc.pipeline import train_rf as train_rf_imported
from gnomad.utils.file_utils import file_exists
from gnomad.resources.resource_utils import TableResource, MatrixTableResource
from gnomad.utils.filtering import add_filters_expr

from gnomad.variant_qc.random_forest import (
    apply_rf_model,
    load_model,
    median_impute_features,
    pretty_print_runs,
    save_model,
)


os.environ['PYSPARK_PYTHON'] = sys.executable

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

project_root = Path(__file__).parent.parent.parent
print(project_root)

s3credentials = os.path.join(
    project_root, "hail_configuration_files/s3_credentials.json")
print(s3credentials)

storage = os.path.join(project_root, "hail_configuration_files/storage.json")

thresholds = os.path.join(
    project_root, "hail_configuration_files/thresholds.json")

with open(f"{s3credentials}", 'r') as f:
    credentials = json.load(f)

with open(f"{storage}", 'r') as f:
    storage = json.load(f)

with open(f"{thresholds}", 'r') as f:
    thresholds = json.load(f)


TRUTH_DATA = ["hapmap", "omni", "mills", "kgp_phase1_hc"]
INBREEDING_COEFF_HARD_CUTOFF = -0.3
tmp_dir = "hdfs://spark-master:9820/"
temp_dir = "file:///home/ubuntu/data/tmp"
#plot_dir = "/home/ubuntu/data/tmp"
plot_dir="/lustre/scratch123/teams/hgi/mercury/megaWES-variantqc"
lustre_dir = "file:///lustre/scratch123/teams/hgi/mercury/megaWES-variantqc"




def add_rank(
    ht: hl.Table,
    score_expr: hl.expr.NumericExpression,
    subrank_expr: Optional[Dict[str, hl.expr.BooleanExpression]] = None,
) -> hl.Table:
    """
    Adds rank based on the `score_expr`. Rank is added for snvs and indels separately.
    If one or more `subrank_expr` are provided, then subrank is added based on all sites for which the boolean expression is true.

    In addition, variant counts (snv, indel separately) is added as a global (`rank_variant_counts`).

    :param ht: input Hail Table containing variants (with QC annotations) to be ranked
    :param score_expr: the Table annotation by which ranking should be scored
    :param subrank_expr: Any subranking to be added in the form name_of_subrank: subrank_filtering_expr
    :return: Table with rankings added
    """

    key = ht.key
    if subrank_expr is None:
        subrank_expr = {}

    temp_expr = {"_score": score_expr}
    temp_expr.update({f"_{name}": expr for name, expr in subrank_expr.items()})
    rank_ht = ht.select(
        **temp_expr, is_snv=hl.is_snp(ht.alleles[0], ht.alleles[1]))

    rank_ht = rank_ht.key_by("_score").persist()
    scan_expr = {
        "rank": hl.cond(
            rank_ht.is_snv,
            hl.scan.count_where(rank_ht.is_snv),
            hl.scan.count_where(~rank_ht.is_snv),
        )
    }
    scan_expr.update(
        {
            name: hl.or_missing(
                rank_ht[f"_{name}"],
                hl.cond(
                    rank_ht.is_snv,
                    hl.scan.count_where(rank_ht.is_snv & rank_ht[f"_{name}"]),
                    hl.scan.count_where(~rank_ht.is_snv & rank_ht[f"_{name}"]),
                ),
            )
            for name in subrank_expr
        }
    )
    rank_ht = rank_ht.annotate(**scan_expr)

    rank_ht = rank_ht.key_by(*key).persist()
    rank_ht = rank_ht.select(*scan_expr.keys())

    ht = ht.annotate(**rank_ht[key])
    return ht


def create_binned_data_initial(ht: hl.Table, data: str, data_type: str, n_bins: int) -> hl.Table:
    # Count variants for ranking
    count_expr = {x: hl.agg.filter(hl.is_defined(ht[x]), hl.agg.counter(hl.cond(hl.is_snp(
        ht.alleles[0], ht.alleles[1]), 'snv', 'indel'))) for x in ht.row if x.endswith('rank')}
    rank_variant_counts = ht.aggregate(hl.Struct(**count_expr))
    logger.info(
        f"Found the following variant counts:\n {pformat(rank_variant_counts)}")
    ht_truth_data = hl.read_table(
        f"{lustre_dir}/variant_qc/truthset_table.ht")
    ht = ht.annotate_globals(rank_variant_counts=rank_variant_counts)
    ht = ht.annotate(
        **ht_truth_data[ht.key],
        # **fam_ht[ht.key],
        # **gnomad_ht[ht.key],
        # **denovo_ht[ht.key],
        # clinvar=hl.is_defined(clinvar_ht[ht.key]),
        indel_length=hl.abs(ht.alleles[0].length()-ht.alleles[1].length()),
        rank_bins=hl.array(
            [hl.Struct(
                rank_id=rank_name,
                bin=hl.int(hl.ceil(hl.float(ht[rank_name] + 1) / hl.floor(ht.globals.rank_variant_counts[rank_name][hl.cond(
                    hl.is_snp(ht.alleles[0], ht.alleles[1]), 'snv', 'indel')] / n_bins)))
            )
                for rank_name in rank_variant_counts]
        ),
        # lcr=hl.is_defined(lcr_intervals[ht.locus])
    )

    ht = ht.explode(ht.rank_bins)
    ht = ht.transmute(
        rank_id=ht.rank_bins.rank_id,
        bin=ht.rank_bins.bin
    )
    ht = ht.filter(hl.is_defined(ht.bin))

    ht = ht.checkpoint(
        f'{lustre_dir}/gnomad_score_binning_tmp.ht', overwrite=True)

    # Create binned data
    return (
        ht
        .group_by(
            rank_id=ht.rank_id,
            contig=ht.locus.contig,
            snv=hl.is_snp(ht.alleles[0], ht.alleles[1]),
            bi_allelic=hl.is_defined(ht.biallelic_rank),
            singleton=ht.transmitted_singleton,
            trans_singletons=hl.is_defined(ht.singleton_rank),
            de_novo_high_quality=ht.de_novo_high_quality_rank,
            de_novo_medium_quality=hl.is_defined(
                ht.de_novo_medium_quality_rank),
            de_novo_synonymous=hl.is_defined(ht.de_novo_synonymous_rank),
            # release_adj=ht.ac > 0,
            bin=ht.bin
        )._set_buffer_size(20000)
        .aggregate(
            min_score=hl.agg.min(ht.score),
            max_score=hl.agg.max(ht.score),
            n=hl.agg.count(),
            n_ins=hl.agg.count_where(
                hl.is_insertion(ht.alleles[0], ht.alleles[1])),
            n_del=hl.agg.count_where(
                hl.is_deletion(ht.alleles[0], ht.alleles[1])),
            n_ti=hl.agg.count_where(hl.is_transition(
                ht.alleles[0], ht.alleles[1])),
            n_tv=hl.agg.count_where(hl.is_transversion(
                ht.alleles[0], ht.alleles[1])),
            n_1bp_indel=hl.agg.count_where(ht.indel_length == 1),
            n_mod3bp_indel=hl.agg.count_where((ht.indel_length % 3) == 0),
            # n_clinvar=hl.agg.count_where(ht.clinvar),
            n_singleton=hl.agg.count_where(ht.transmitted_singleton),
            n_high_quality_de_novos=hl.agg.count_where(
                ht.de_novo_data.p_de_novo[0] > 0.99),
            #n_validated_DDD_denovos=hl.agg.count_where(
            #    ht.inheritance.contains("De novo")),
            n_medium_quality_de_novos=hl.agg.count_where(
                ht.de_novo_data.p_de_novo[0] > 0.5),
            n_high_confidence_de_novos=hl.agg.count_where(
                ht.de_novo_data.confidence[0] == 'HIGH'),
            n_de_novo=hl.agg.filter(ht.family_stats.unrelated_qc_callstats.AC[0][1] == 0, hl.agg.sum(
                ht.family_stats.mendel[0].errors)),
            n_high_quality_de_novos_synonymous=hl.agg.count_where(
                (ht.de_novo_data.p_de_novo[0] > 0.99) & (ht.consequence == "synonymous_variant")),
            n_trans_singletons_synonymous_algorithm=hl.agg.count_where(
                ht.variant_transmitted_singletons > 0 ),
            n_untrans_singletons_synonymous_algorithm=hl.agg.count_where(
                ht.variant_untransmitted_singletons > 0),
            validated_de_novos=hl.agg.count_where(ht.validated_denovo_inheritance=="De novo constitutive"),
            # n_de_novo_no_lcr=hl.agg.filter(~ht.lcr & (
            #    ht.family_stats.unrelated_qc_callstats.AC[1] == 0), hl.agg.sum(ht.family_stats.mendel.errors)),
            n_de_novo_sites=hl.agg.filter(ht.family_stats.unrelated_qc_callstats.AC[0][1] == 0, hl.agg.count_where(
                ht.family_stats.mendel[0].errors > 0)),
            # n_de_novo_sites_no_lcr=hl.agg.filter(~ht.lcr & (
            #    ht.family_stats.unrelated_qc_callstats.AC[1] == 0), hl.agg.count_where(ht.family_stats.mendel.errors > 0)),
            n_trans_singletons=hl.agg.filter((ht.ac_raw < 3) & (
                ht.family_stats.unrelated_qc_callstats.AC[0][1] == 1), hl.agg.sum(ht.family_stats.tdt[0].t)),
            n_trans_singletons_synonymous_hail=hl.agg.filter((ht.ac_raw < 3) & (ht.consequence == "synonymous_variant") & (
                ht.family_stats.unrelated_qc_callstats.AC[0][1] == 1), hl.agg.sum(ht.family_stats.tdt[0].t)),
            n_untrans_singletons=hl.agg.filter((ht.ac_raw < 3) & (
                ht.family_stats.unrelated_qc_callstats.AC[0][1] == 1), hl.agg.sum(ht.family_stats.tdt[0].u)),
            n_untrans_singletons_synonymous_hail=hl.agg.filter((ht.ac_raw < 3) & (ht.consequence == "synonymous_variant") & (
                ht.family_stats.unrelated_qc_callstats.AC[0][1] == 1), hl.agg.sum(ht.family_stats.tdt[0].u)),
            n_train_trans_singletons=hl.agg.count_where(
                (ht.family_stats.unrelated_qc_callstats.AC[0][1] == 1) & (ht.family_stats.tdt[0].t == 1)),
            n_omni=hl.agg.count_where(ht.omni),
            n_mills=hl.agg.count_where(ht.mills),
            n_hapmap=hl.agg.count_where(ht.hapmap),
            n_kgp_high_conf_snvs=hl.agg.count_where(
                ht.kgp_phase1_hc),
            fail_hard_filters=hl.agg.count_where(ht.fail_hard_filters),
            # n_vqsr_pos_train=hl.agg.count_where(ht.vqsr_positive_train_site),
            # n_vqsr_neg_train=hl.agg.count_where(ht.vqsr_negative_train_site)
        )
    )


######################################
# main
########################################


def main(args):

    # ht after random model
    run_hash = args.run_hash
    ht = hl.read_table(f'{lustre_dir}/variant_qc/models/{run_hash}_rf_result_FINAL_for_RANKING_denovo_validated.ht')

    if args.add_rank:
        ht_ranked = add_rank(ht,
                             score_expr=1-ht.rf_probability["FP"],
                             #score_expr=ht.rf_probability["TP"],
                             subrank_expr={
                                 'singleton_rank': ht.transmitted_singleton,
                                 'biallelic_rank': ~ht.was_split,
                                 'biallelic_singleton_rank': ~ht.was_split & ht.transmitted_singleton,
                                 'de_novo_high_quality_rank': ht.de_novo_data.p_de_novo[0] > 0.9,
                                 'de_novo_medium_quality_rank': ht.de_novo_data.p_de_novo[0] > 0.5,
                                 'de_novo_synonymous_rank': ht.consequence == "synonymous_variant",
                                 # 'adj_rank': ht.ac > 0,
                                 # 'adj_biallelic_rank': ~ht.was_split & (ht.ac > 0),
                                 # 'adj_singleton_rank': ht.transmitted_singleton & (ht.ac > 0),
                                 # 'adj_biallelic_singleton_rank': ~ht.was_split & ht.transmitted_singleton & (ht.ac > 0)
                             }
                             )
        ht_ranked = ht_ranked.annotate(score=1-ht_ranked.rf_probability["TP"])
        #ht_ranked = ht_ranked.annotate(score=ht_ranked.rf_probability["TP"])
        ht_ranked = ht_ranked.checkpoint(
            f'{lustre_dir}/variant_qc/models/{run_hash}_rf_result_ranked.ht', overwrite=True)

    if args.add_bin:

        # ht = hl.read_table(
        #    f'{lustre_dir}/ddd-elgh-ukbb/variant_qc/models/{run_hash}/{run_hash}_rf_result_ranked.ht')
        ht = hl.read_table(
            f'{lustre_dir}/variant_qc/models/{run_hash}_rf_result_ranked.ht')

       
        ht_bins = create_binned_data_initial(ht, "exomes", "RF", n_bins=100)
        ht_bins.write(
            f'{lustre_dir}/variant_qc/models/{run_hash}_rf_result_ranked_BINS.ht', overwrite=True)
        # ht_grouped = compute_grouped_binned_ht(ht_bins)
        # ht_grouped.write(
        #    f'{tmp_dir}/ddd-elgh-ukbb/{run_hash}_rf_result_ranked_BINS_Grouped.ht', overwrite=True)


if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory

    hl.init(sc=sc, tmp_dir=lustre_dir, local_tmpdir=lustre_dir, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])
    n_partitions = 500
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--run_hash",
        help="Run hash. Created by --train_rf and only needed for --apply_rf without running --train_rf",
        required=False,
    )

    actions = parser.add_argument_group("Actions")

    actions.add_argument(
        "--add_rank",
        help="Add rank to RF results",
        action="store_true",
    )
    actions.add_argument(
        "--add_bin",
        help="Split to bin and calculate stats for RF results",
        action="store_true",
    )
    rf_params = parser.add_argument_group("Random Forest Parameters")
    rf_params.add_argument(
        "--fp_to_tp",
        help="Ratio of FPs to TPs for training the RF model. If 0, all training examples are used. (default=1.0)",
        default=1.0,
        type=float,
    )
    rf_params.add_argument(
        "--test_intervals",
        help='The specified interval(s) will be held out for testing and evaluation only. (default to "chr20")',
        nargs="+",
        type=str,
        default="chr4",
    )
    rf_params.add_argument(
        "--num_trees",
        help="Number of trees in the RF model. (default=500)",
        default=500,
        type=int,
    )
    rf_params.add_argument(
        "--max_depth",
        help="Maxmimum tree depth in the RF model. (default=5)",
        default=5,
        type=int,
    )
    training_params = parser.add_argument_group("Training data parameters")
    training_params.add_argument(
        "--adj", help="Use adj genotypes.", action="store_true"
    )
    training_params.add_argument(
        "--vqsr_training", help="Use VQSR training examples", action="store_true"
    )
    training_params.add_argument(
        "--vqsr_type",
        help="If a string is provided the VQSR training annotations will be used for training.",
        default="alleleSpecificTrans",
        choices=["classic", "alleleSpecific", "alleleSpecificTrans"],
        type=str,
    )
    training_params.add_argument(
        "--no_transmitted_singletons",
        help="Do not use transmitted singletons for training.",
        action="store_true",
    )
    training_params.add_argument(
        "--no_inbreeding_coeff",
        help="Train RF without inbreeding coefficient as a feature.",
        action="store_true",
    )

    finalize_params = parser.add_argument_group("Finalize RF Table parameters")
    finalize_params.add_argument(
        "--snp_cutoff", help="Percentile to set RF cutoff", type=float, default=90.0
    )
    finalize_params.add_argument(
        "--indel_cutoff", help="Percentile to set RF cutoff", type=float, default=80.0
    )
    finalize_params.add_argument(
        "--treat_cutoff_as_prob",
        help="If set snp_cutoff and indel_cutoff will be probability rather than percentile ",
        action="store_true",
    )
    args = parser.parse_args()
    main(args)
