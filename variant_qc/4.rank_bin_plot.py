# Pavlos Antoniou
# 16/09/2020
#  trio matrixtable creation from fam file
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
#from gnomad.variant_qc.pipeline import train_rf as train_rf_imported
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
        "rank": hl.if_else(
            rank_ht.is_snv,
            hl.scan.count_where(rank_ht.is_snv),
            hl.scan.count_where(~rank_ht.is_snv),
        )
    }
    scan_expr.update(
        {
            name: hl.or_missing(
                rank_ht[f"_{name}"],
                hl.if_else(
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


def compute_quantile_bin(
    ht: hl.Table,
    score_expr: hl.expr.NumericExpression,
    bin_expr: Dict[str, hl.expr.BooleanExpression] = {"bin": True},
    compute_snv_indel_separately: bool = True,
    n_bins: int = 100,
    k: int = 1000,
    desc: bool = True,
) -> hl.Table:
    """
    Returns a table with a bin for each row based on quantiles of `score_expr`.
    The bin is computed by dividing the `score_expr` into `n_bins` bins containing an equal number of elements.
    This is done based on quantiles computed with hl.agg.approx_quantiles. If a single value in `score_expr` spans more
    than one bin, the rows with this value are distributed randomly across the bins it spans.
    If `compute_snv_indel_separately` is True all items in `bin_expr` will be stratified by snv / indels for the bin
    calculation. Because SNV and indel rows are mutually exclusive, they are re-combined into a single annotation. For
    example if we have the following four variants and scores and `n_bins` of 2:
    ========   =======   ======   =================   =================
    Variant    Type      Score    bin - `compute_snv_indel_separately`:
    --------   -------   ------   -------------------------------------
    \          \         \        False               True
    ========   =======   ======   =================   =================
    Var1       SNV       0.1      1                   1
    Var2       SNV       0.2      1                   2
    Var3       Indel     0.3      2                   1
    Var4       Indel     0.4      2                   2
    ========   =======   ======   =================   =================
    .. note::
        The `bin_expr` defines which data the bin(s) should be computed on. E.g., to get a biallelic quantile bin and an
        singleton quantile bin, the following could be used:
        .. code-block:: python
            bin_expr={
                'biallelic_bin': ~ht.was_split,
                'singleton_bin': ht.singleton
            }
    :param ht: Input Table
    :param score_expr: Expression containing the score
    :param bin_expr: Quantile bin(s) to be computed (see notes)
    :param compute_snv_indel_separately: Should all `bin_expr` items be stratified by snv / indels
    :param n_bins: Number of bins to bin the data into
    :param k: The `k` parameter of approx_quantiles
    :param desc: Whether to bin the score in descending order
    :return: Table with the quantile bins
    """
    import math

    def quantiles_to_bin_boundaries(quantiles: List[int]) -> Dict:
        """
        Merges bins with the same boundaries into a unique bin while keeping track of
        which bins have been merged and the global index of all bins.
        :param quantiles: Original bins boundaries
        :return: (dict of the indices of bins for which multiple bins were collapsed -> number of bins collapsed,
                  Global indices of merged bins,
                  Merged bins boundaries)
        """

        # Pad the quantiles to create boundaries for the first and last bins
        bin_boundaries = [-math.inf] + quantiles + [math.inf]
        merged_bins = defaultdict(int)

        # If every quantile has a unique value, then bin boudaries are unique
        # and can be passed to binary_search as-is
        if len(quantiles) == len(set(quantiles)):
            return dict(
                merged_bins=merged_bins,
                global_bin_indices=list(range(len(bin_boundaries))),
                bin_boundaries=bin_boundaries,
            )

        indexed_bins = list(enumerate(bin_boundaries))
        i = 1
        while i < len(indexed_bins):
            if indexed_bins[i - 1][1] == indexed_bins[i][1]:
                merged_bins[i - 1] += 1
                indexed_bins.pop(i)
            else:
                i += 1

        return dict(
            merged_bins=merged_bins,
            global_bin_indices=[x[0] for x in indexed_bins],
            bin_boundaries=[x[1] for x in indexed_bins],
        )

    if compute_snv_indel_separately:
        # For each bin, add a SNV / indel stratification
        bin_expr = {
            f"{bin_id}_{snv}": (bin_expr & snv_expr)
            for bin_id, bin_expr in bin_expr.items()
            for snv, snv_expr in [
                ("snv", hl.is_snp(ht.alleles[0], ht.alleles[1])),
                ("indel", ~hl.is_snp(ht.alleles[0], ht.alleles[1])),
            ]
        }
        print("ADSADSADASDAS")
        print(bin_expr)

    bin_ht = ht.annotate(
        **{f"_filter_{bin_id}": bin_expr for bin_id, bin_expr in bin_expr.items()},
        _score=score_expr,
        snv=hl.is_snp(ht.alleles[0], ht.alleles[1]),
    )
    print(bin_ht.show())
    logger.info(
        f"Adding quantile bins using approximate_quantiles binned into {n_bins}, using k={k}"
    )
    bin_stats = bin_ht.aggregate(
        hl.struct(
            **{
                bin_id: hl.agg.filter(
                    bin_ht[f"_filter_{bin_id}"],
                    hl.struct(
                        n=hl.agg.count(),
                        quantiles=hl.agg.approx_quantiles(
                            bin_ht._score, [x / (n_bins) for x in range(1, n_bins)], k=k
                        ),
                    ),
                )
                for bin_id in bin_expr
            }
        )
    )

    # Take care of bins with duplicated boundaries
    bin_stats = bin_stats.annotate(
        **{
            rname: bin_stats[rname].annotate(
                **quantiles_to_bin_boundaries(bin_stats[rname].quantiles)
            )
            for rname in bin_stats
        }
    )

    bin_ht = bin_ht.annotate_globals(
        bin_stats=hl.literal(
            bin_stats,
            dtype=hl.tstruct(
                **{
                    bin_id: hl.tstruct(
                        n=hl.tint64,
                        quantiles=hl.tarray(hl.tfloat64),
                        bin_boundaries=hl.tarray(hl.tfloat64),
                        global_bin_indices=hl.tarray(hl.tint32),
                        merged_bins=hl.tdict(hl.tint32, hl.tint32),
                    )
                    for bin_id in bin_expr
                }
            ),
        )
    )

    # Annotate the bin as the index in the unique boundaries array
    bin_ht = bin_ht.annotate(
        **{
            bin_id: hl.or_missing(
                bin_ht[f"_filter_{bin_id}"],
                hl.binary_search(
                    bin_ht.bin_stats[bin_id].bin_boundaries, bin_ht._score
                ),
            )
            for bin_id in bin_expr
        }
    )

    # Convert the bin to global bin by expanding merged bins, that is:
    # If a value falls in a bin that needs expansion, assign it randomly to one of the expanded bins
    # Otherwise, simply modify the bin to its global index (with expanded bins that is)
    bin_ht = bin_ht.select(
        "snv",
        **{
            bin_id: hl.if_else(
                bin_ht.bin_stats[bin_id].merged_bins.contains(bin_ht[bin_id]),
                bin_ht.bin_stats[bin_id].global_bin_indices[bin_ht[bin_id]]
                + hl.int(
                    hl.rand_unif(
                        0, bin_ht.bin_stats[bin_id].merged_bins[bin_ht[bin_id]] + 1
                    )
                ),
                bin_ht.bin_stats[bin_id].global_bin_indices[bin_ht[bin_id]],
            )
            for bin_id in bin_expr
        },
    )

    if desc:
        bin_ht = bin_ht.annotate(
            **{bin_id: n_bins - bin_ht[bin_id] for bin_id in bin_expr}
        )

    # Because SNV and indel rows are mutually exclusive, re-combine them into a single bin.
    # Update the global bin_stats struct to reflect the change in bin names in the table
    if compute_snv_indel_separately:
        bin_expr_no_snv = {bin_id.rsplit(
            "_", 1)[0] for bin_id in bin_ht.bin_stats}
        bin_ht = bin_ht.annotate_globals(
            bin_stats=hl.struct(
                **{
                    bin_id: hl.struct(
                        **{
                            snv: bin_ht.bin_stats[f"{bin_id}_{snv}"]
                            for snv in ["snv", "indel"]
                        }
                    )
                    for bin_id in bin_expr_no_snv
                }
            )
        )

        bin_ht = bin_ht.transmute(
            **{
                bin_id: hl.if_else(
                    bin_ht.snv, bin_ht[f"{bin_id}_snv"], bin_ht[f"{bin_id}_indel"],
                )
                for bin_id in bin_expr_no_snv
            }
        )

    return bin_ht

######################################
# main
########################################


def main(args):

    # ht after random model
    run_hash = args.run_hash
    ht = hl.read_table(
        f'{temp_dir}/ddd-elgh-ukbb/variant_qc/models/{run_hash}/rf_result_sanger_cohorts_new.ht')

    if args.add_rank:
        ht_ranked = add_rank(ht, 'rf_probability')
        ht_ranked = ht_ranked.checkpoint(
            f'{tmp_dir}/ddd-elgh-ukbb/{run_hash}_rf_result_ranked.ht', overwrite=True)

    if args.add_bin:
        ht = hl.read_table(
            f'{temp_dir}/ddd-elgh-ukbb/variant_qc/models/{run_hash}/{run_hash}_rf_result_ranked.ht')
        ht_bins = compute_quantile_bin(ht, ht.InbreedingCoeff, bin_expr={
            'biallelic_bin': ~ht.was_split,
            'singleton_bin': ht.transmitted_singleton,
        }, compute_snv_indel_separately=True, n_bins=100, k=500, desc=True)
        ht_bins.write(
            f'{tmp_dir}/ddd-elgh-ukbb/{run_hash}_rf_result_ranked_BINS.ht', overwrite=True)


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
