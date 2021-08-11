# Pavlos Antoniou
# 16/09/2020
#  trio matrixtable creation from fam file
from hail import Table
import os
import argparse
import hail as hl
import pandas as pd
import numpy as np
import pyspark
from pprint import pformat
import json
import sys
import re
from pathlib import Path
import logging
from typing import Any, Counter, List, Optional, Tuple, Union, Dict
import uuid
import json
from bokeh.plotting import output_file, save, show
from gnomad.resources.grch38 import gnomad
from gnomad.utils.annotations import unphase_call_expr, add_variant_type
from gnomad.variant_qc.pipeline import create_binned_ht, score_bin_agg, train_rf_model
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

project_root = Path(__file__).parent.parent
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

LABEL_COL = "rf_label"
TRAIN_COL = "rf_train"
PREDICTION_COL = "rf_prediction"
INFO_FEATURES = [
    "AS_QD",
    "AS_ReadPosRankSum",
    "AS_MQRankSum",
    "AS_SOR",
]  # Note: AS_SOR is currently in VQSR HT and named SOR in the VQSR split HT
FEATURES = [
    "InbreedingCoeff",
    "variant_type",
    "allele_type",
    "n_alt_alleles",
    "was_mixed",
    "has_star",
    # "AS_QD",
    #    "AS_MQRankSum",
    #    "AS_SOR",
    #    "AS_ReadPosRankSum",
]

TRUTH_DATA = ["hapmap", "omni", "mills", "kgp_phase1_hc"]
INBREEDING_COEFF_HARD_CUTOFF = -0.3
tmp_dir = "hdfs://spark-master:9820/"
temp_dir = "file:///home/ubuntu/data/tmp"
plot_dir = "/home/ubuntu/data/tmp"


def get_rf(
    data: str = "rf_result",
    run_hash: Optional[str] = None,
) -> Union[str, TableResource]:
    """
    Gets the path to the desired RF data.
    Data can take the following values:
        - 'training': path to the training data for a given run
        - 'model': path to pyspark pipeline RF model
        - 'rf_result' (default): path to HT containing result of RF filtering
    :param str data: One of 'training', 'model' or 'rf_result' (default)
    :param str run_hash: Hash of RF run to load
    :return: Path to desired RF data
    """

    if data == "model":
        return f"{tmp_dir}/models/{run_hash}/{data}.model"
    else:
        return TableResource(f"{tmp_dir}/models/{run_hash}/{data}.ht")


def get_rf_runs(rf_json_fp: str) -> Dict:
    """
    Loads RF run data from JSON file.
    :param rf_json_fp: File path to rf json file.
    :return: Dictionary containing the content of the JSON file, or an empty dictionary if the file wasn't found.
    """
    if file_exists(rf_json_fp):
        with hl.hadoop_open(rf_json_fp) as f:
            return json.load(f)
    else:
        logger.warning(
            f"File {rf_json_fp} could not be found. Returning empty RF run hash dict."
        )
        return {}


def get_truth_ht() -> Table:
    """
    Returns a table with the following annotations from the latest version of the corresponding truth data:
    - hapmap
    - kgp_omni (1000 Genomes intersection Onni 2.5M array)
    - kgp_phase_1_hc (high confidence sites in 1000 genonmes)
    - mills (Mills & Devine indels)
    :return: A table with the latest version of popular truth data annotations
    """
    omni_ht = hl.read_table(omni)
    mills_ht = hl.read_table(mills)
    thousand_genomes_ht = hl.read_table(thousand_genomes)
    hapmap_ht = hl.read_table(hapmap)
    return (
        hapmap_ht
        .select(hapmap=True)
        .join(omni_ht.select(omni=True), how="outer")
        .join(thousand_genomes_ht.select(kgp_phase1_hc=True), how="outer")
        .join(mills_ht.select(mills=True), how="outer")
        .repartition(200, shuffle=False)
        .persist()
    )


def generate_allele_data(mt: hl.MatrixTable) -> hl.Table:
    """
    Writes bi-allelic sites MT with the following annotations:
     - allele_data (nonsplit_alleles, has_star, variant_type, and n_alt_alleles)

    :param MatrixTable mt: Full unsplit MT
    :return: Table with allele data annotations
    :rtype: Table
    """
    ht = mt.rows().select()
    allele_data = hl.struct(nonsplit_alleles=ht.alleles,
                            has_star=hl.any(lambda a: a == '*', ht.alleles))
    ht = ht.annotate(allele_data=allele_data.annotate(
        **add_variant_type(ht.alleles)))

    ht = hl.split_multi_hts(ht)
    allele_type = (hl.case()
                   .when(hl.is_snp(ht.alleles[0], ht.alleles[1]), 'snv')
                   .when(hl.is_insertion(ht.alleles[0], ht.alleles[1]), 'ins')
                   .when(hl.is_deletion(ht.alleles[0], ht.alleles[1]), 'del')
                   .default('complex')
                   )
    ht = ht.annotate(allele_data=ht.allele_data.annotate(allele_type=allele_type,
                                                         was_mixed=ht.allele_data.variant_type == 'mixed'))
    return ht


def save_model(
    rf_pipeline: pyspark.ml.PipelineModel, out_path: str, overwrite: bool = False
) -> None:
    """
    Saves a Random Forest pipeline model.
    :param rf_pipeline: Pipeline to save
    :param out_path: Output path
    :param overwrite: If set, will overwrite existing file(s) at output location
    :return: Nothing
    """
    logger.info("Saving model to %s" % out_path)
    if overwrite:
        rf_pipeline.write().overwrite().save(out_path)
    else:
        rf_pipeline.save(out_path)


def create_grouped_bin_ht(model_id: str, overwrite: bool = False) -> None:
    """
    Creates binned data from a quantile bin annotated Table grouped by bin_id (rank, bi-allelic, etc.), contig, snv,
    bi_allelic and singleton containing the information needed for evaluation plots.
    :param str model_id: Which data/run hash is being created
    :param bool overwrite: Should output files be overwritten if present
    :return: None
    :rtype: None
    """

    ht = get_score_quantile_bins(model_id, aggregated=False).ht()

    # Count variants for ranking
    count_expr = {
        x: hl.agg.filter(
            hl.is_defined(ht[x]),
            hl.agg.counter(
                hl.cond(hl.is_snp(ht.alleles[0],
                                  ht.alleles[1]), "snv", "indel")
            ),
        )
        for x in ht.row
        if x.endswith("bin")
    }
    bin_variant_counts = ht.aggregate(hl.struct(**count_expr))
    logger.info(
        f"Found the following variant counts:\n {pformat(bin_variant_counts)}")
    ht = ht.annotate_globals(bin_variant_counts=bin_variant_counts)

    trio_stats_ht = fam_stats.ht()

    logger.info(f"Creating grouped bin table...")
    grouped_binned_ht = compute_grouped_binned_ht(
        ht, checkpoint_path=get_checkpoint_path(f"grouped_bin_{model_id}"),
    )

    logger.info(f"Aggregating grouped bin table...")
    agg_ht = grouped_binned_ht.aggregate(
        **score_bin_agg(grouped_binned_ht, fam_stats_ht=trio_stats_ht)
    )

    agg_ht.write(
        get_score_quantile_bins(
            model_id, aggregated=True).path, overwrite=overwrite,
    )


def train_rf(ht, args):
    features = FEATURES
    test_intervals = args.test_intervals
    if args.no_inbreeding_coeff:
        features.remove("InbreedingCoeff")

    fp_expr = ht.fail_hard_filters
    tp_expr = ht.omni | ht.mills
    if not args.no_transmitted_singletons:
        tp_expr = tp_expr | ht.transmitted_singleton

    if test_intervals:
        if isinstance(test_intervals, str):
            test_intervals = [test_intervals]
        test_intervals = [
            hl.parse_locus_interval(x, reference_genome="GRCh38")
            for x in test_intervals
        ]

    ht = ht.annotate(tp=tp_expr, fp=fp_expr)

    rf_ht, rf_model = train_rf_model(
        ht,
        rf_features=features,
        tp_expr=ht.tp,
        fp_expr=ht.fp,
        fp_to_tp=args.fp_to_tp,
        num_trees=args.num_trees,
        max_depth=args.max_depth,
        test_expr=hl.literal(test_intervals).any(
            lambda interval: interval.contains(ht.locus)
        ),
    )

    logger.info("Joining original RF Table with training information")
    ht = ht.join(rf_ht, how="left")

    return ht, rf_model


def get_run_data(
    transmitted_singletons: bool,
    adj: bool,
    vqsr_training: bool,
    test_intervals: List[str],
    features_importance: Dict[str, float],
    test_results: List[hl.tstruct],
) -> Dict:
    """
    Creates a Dict containing information about the RF input arguments and feature importance
    :param bool transmitted_singletons: True if transmitted singletons were used in training
    :param bool adj: True if training variants were filtered by adj
    :param bool vqsr_training: True if VQSR training examples were used for RF training
    :param List of str test_intervals: Intervals withheld from training to be used in testing
    :param Dict of float keyed by str features_importance: Feature importance returned by the RF
    :param List of struct test_results: Accuracy results from applying RF model to the test intervals
    :return: Dict of RF information
    """
    if vqsr_training:
        transmitted_singletons = None

    run_data = {
        "input_args": {
            "transmitted_singletons": transmitted_singletons,
            "adj": adj,
            "vqsr_training": vqsr_training,
        },
        "features_importance": features_importance,
        "test_intervals": test_intervals,
    }

    if test_results is not None:
        tps = 0
        total = 0
        for row in test_results:
            values = list(row.values())
            # Note: values[0] is the TP/FP label and values[1] is the prediction
            if values[0] == values[1]:
                tps += values[2]
            total += values[2]
        run_data["test_results"] = [dict(x) for x in test_results]
        run_data["test_accuracy"] = tps / total

    return run_data


def get_score_quantile_bins(model_id: str, aggregated: bool) -> TableResource:
    return TableResource('{}/{}.{}.ht'.format(
        f"{tmp_dir}",
        model_id,
        'binned' if aggregated else 'rank'
    ))


def generate_final_rf_ht(
    ht: hl.Table,
    snp_cutoff: Union[int, float],
    indel_cutoff: Union[int, float],
    inbreeding_coeff_cutoff: float = INBREEDING_COEFF_HARD_CUTOFF,
    determine_cutoff_from_bin: bool = False,
    aggregated_bin_ht: Optional[hl.Table] = None,
    bin_id: Optional[hl.expr.Int32Expression] = None,
) -> hl.Table:
    """
    Prepares finalized RF model given an RF result table from `rf.apply_rf_model` and cutoffs for filtering.
    If `determine_cutoff_from_bin` is True, `aggregated_bin_ht` must be supplied to determine the SNP and indel RF
    probabilities to use as cutoffs from an aggregated quantile bin Table like one created by
    `compute_grouped_binned_ht` in combination with `score_bin_agg`.
    :param ht: RF result table from `rf.apply_rf_model` to prepare as the final RF Table
    :param ac0_filter_expr: Expression that indicates if a variant should be filtered as allele count 0 (AC0)
    :param ts_ac_filter_expr: Expression in `ht` that indicates if a variant is a transmitted singleton
    :param mono_allelic_fiter_expr: Expression indicating if a variant is mono-allelic
    :param snp_cutoff: RF probability or bin (if `determine_cutoff_from_bin` True) to use for SNP variant QC filter
    :param indel_cutoff: RF probability or bin (if `determine_cutoff_from_bin` True) to use for indel variant QC filter
    :param inbreeding_coeff_cutoff: InbreedingCoeff hard filter to use for variants
    :param determine_cutoff_from_bin: If True RF probability will be determined using bin info in `aggregated_bin_ht`
    :param aggregated_bin_ht: File with aggregate counts of variants based on quantile bins
    :param bin_id: Name of bin to use in 'bin_id' column of `aggregated_bin_ht` to use to determine probability cutoff
    :return: Finalized random forest Table annotated with variant filters
    """
    # Determine SNP and indel RF cutoffs if given bin instead of RF probability

    snp_cutoff_global = hl.struct(min_score=snp_cutoff)
    indel_cutoff_global = hl.struct(min_score=indel_cutoff)

    # Add filters to RF HT
    filters = dict()

    if ht.any(hl.is_missing(ht.rf_probability["TP"])):
        raise ValueError("Missing RF probability!")

    filters["RF"] = (
        hl.is_snp(ht.alleles[0], ht.alleles[1])
        & (ht.rf_probability["TP"] < snp_cutoff_global.min_score)
    ) | (
        ~hl.is_snp(ht.alleles[0], ht.alleles[1])
        & (ht.rf_probability["TP"] < indel_cutoff_global.min_score)
    )

    # Fix annotations for release
    annotations_expr = {
        "rf_positive_label": hl.or_else(ht.tp, False),
        "rf_negative_label": ht.fail_hard_filters,
        "rf_probability": ht.rf_probability["TP"], }

    ht = ht.transmute(filters=add_filters_expr(
        filters=filters), **annotations_expr)

    ht = ht.annotate_globals(
        rf_snv_cutoff=snp_cutoff_global, rf_indel_cutoff=indel_cutoff_global
    )

    return ht


######################################
# main
########################################

# annotate with RF probability score 1-
    # 'rf_probability': dict<str, float64>
    #  'rf_train': bool
    # 'rf_label': str
    # 'rf_test': bool
    # 'rf_probability': dict<str, float64> :  - rf_probability['TP']
    # 'rf_prediction': str
    # 'consequence': str
    # 'inheritance': str
    # 'rank': int64
    # 'score': float64

    # info rf_probability['TP']
    # FILTER column to PASS for ≤80 InDels ≤90 SNVs
   # ht = generate_final_rf_ht(
   #     ht,
   #     snp_cutoff=args.snp_cutoff,
   #     indel_cutoff=args.indel_cutoff,
   #     determine_cutoff_from_bin=False,
   #     bin_id=ht.bin,
   #     inbreeding_coeff_cutoff=INBREEDING_COEFF_HARD_CUTOFF,
   # )
    # This column is added by the RF module based on a 0.5 threshold which doesn't correspond to what we use
    #ht = ht.drop(ht[PREDICTION_COL])
def main(args):

    print("main")

    run_hash = "91b132aa"
    ht = hl.read_table(
        f'{temp_dir}/ddd-elgh-ukbb/variant_qc/models/{run_hash}/{run_hash}_rf_result_ranked_denovo_ddd_comp.ht')

    mt = hl.read_matrix_table(
        f'{temp_dir}/ddd-elgh-ukbb/Sanger_cohorts_chr1-7and20_split_sampleqc_filtered.mt')
    mt = mt.annotate_rows(
        Variant_Type=hl.cond((hl.is_snp(mt.alleles[0], mt.alleles[1])), "SNP",
                             hl.cond(
            hl.is_insertion(
                mt.alleles[0], mt.alleles[1]),
            "INDEL",
            hl.cond(hl.is_deletion(mt.alleles[0],
                                   mt.alleles[1]), "INDEL",
                    "Other"))))
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            rf_probability=ht[mt.row_key].rf_probability['TP'])
    )
    mt = mt.annotate_rows(
        info=mt.info.annotate(score=ht[mt.row_key].score)
    )

    filter_column_annotation = (
        hl.case()
        .when(((mt.Variant_Type == "SNP") & (mt.info.rf_probability <= 0.90)), "PASS")
        .when(((mt.Variant_Type == "INDEL") & (mt.info.rf_probability <= 0.80)), "PASS")
        .default(".")  # remove everything else
    )

# mt_annotated = mt.annotate_rows(mt.filters=filter_column_annotation)
    mt1 = mt.annotate_rows(
        filtercol=((filter_column_annotation))
    )
    mt_fail = mt1.filter_rows(mt1.filtercol == ".")
    print(mt_fail.count())

    mt2 = mt1.annotate_rows(filters=mt1.filters.add(mt1.filtercol))
    mt_fail2 = mt2.filter_rows(mt2.filters.contains("."))
    mt_pass = mt2.filter_rows(mt2.filters.contains("PASS"))
    print(mt_fail2.count())
    print(mt_pass.count())

    mt2 = mt2.checkpoint(
        f'{tmp_dir}/Sanger_cohorts_chr1-7and20_after_RF_final.mt', overwrite=True)

    hl.export_vcf(
        mt2, f'{tmp_dir}/Sanger_cohorts_chr1-7and20_after_RF_final.vcf.bgz',parallel='separate_header')


if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory

    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])
    n_partitions = 500
    parser = argparse.ArgumentParser()

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
