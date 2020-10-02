# Pavlos Antoniou
# 16/09/2020
#  trio matrixtable creation from fam file
import os
import argparse
import hail as hl
import pandas as pd
import numpy as np
import pyspark
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
from gnomad.variant_qc.random_forest import (
    apply_rf_model,
    load_model,
    median_impute_features,
    pretty_print_runs,
    save_model,
)


from hail import Table
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


def main(args):

    print("main")
    ht = hl.read_table(
        f'{temp_dir}/ddd-elgh-ukbb/variant_qc/Sanger_table_for_RF_by_variant_type.ht')

    run_hash = str(uuid.uuid4())[:8]
    rf_runs = get_rf_runs(f'{tmp_dir}/ddd-elgh-ukbb/')
    while run_hash in rf_runs:
        run_hash = str(uuid.uuid4())[:8]
    ht_result, rf_model = train_rf(ht, args)
    print("Writing out ht_training data")
    ht_result = ht_result.checkpoint(
        f'{tmp_dir}/ddd-elgh-ukbb/Sanger_RF_training_data.ht', overwrite=True)
    rf_runs[run_hash] = get_run_data(
        transmitted_singletons=True,

        # features_importance=hl.eval(ht.features_importance),
        test_results=hl.eval(ht.test_results),
    )

    with hl.hadoop_open(f'{plot_dir}/ddd-elgh-ukbb/variant_qc/rf_runs.json', "w") as f:
        json.dump(rf_runs, f)

    logger.info("Saving RF model")
    save_model(
        rf_model, f'{tmp_dir}/ddd-elgh-ukbb/rf_model.model')


if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory
    tmp_dir = "hdfs://spark-master:9820/"
    temp_dir = "file:///home/ubuntu/data/tmp"
    plot_dir = "/home/ubuntu/data/tmp"
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])
    n_partitions = 500
    parser = argparse.ArgumentParser()
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
        default="chr20",
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
