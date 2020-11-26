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
from gnomad.variant_qc.pipeline import create_binned_ht, score_bin_agg, train_rf_model
from gnomad.variant_qc.pipeline import test_model, sample_training_examples, get_features_importance
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
# FEATURES = [
#    "InbreedingCoeff",
#    "variant_type",
#    "n_alt_alleles",
#    "allele_type",
#    "has_star",
# "was_mixed",
#    "AS_QD",
# "AS_MQRankSum",
#    "AS_SOR",
#    "AS_ReadPosRankSum",
# ]
FEATURES = [
    "InbreedingCoeff",
    "variant_type",
    "allele_type",
    "n_alt_alleles",
    "was_mixed",
    "has_star",
    "QD",  # QualByDepth
    "MQRankSum",  # MappingQualityRankSumTest
    "SOR",  # StrandOddsRatio
    "ReadPosRankSum",  # ReadPosRankSumTest
    "FS",  # FisherStrand
    "DP"  # Depth
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


def train_rf(ht, args):
    features = FEATURES
    test_intervals = args.test_intervals
    print("test_intervals")
    print(test_intervals)
    #test_intervals = False

    if args.no_inbreeding_coeff:
        features.remove("InbreedingCoeff")

    fp_expr = ht.fail_hard_filters
    tp_expr = ht.omni | ht.mills | ht.kgp_phase1_hc | ht.hapmap
    if not args.no_transmitted_singletons:
        tp_expr = tp_expr | ht.transmitted_singleton

    if test_intervals:

        if isinstance(test_intervals, str):
            test_intervals = [test_intervals]
        test_intervals = [
            hl.parse_locus_interval(x, reference_genome="GRCh38")
            for x in test_intervals
        ]
        print(hl.eval(test_intervals))

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
            lambda interval: interval.contains(ht.locus)),
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
    ac0_filter_expr: hl.expr.BooleanExpression,
    ts_ac_filter_expr: hl.expr.BooleanExpression,
    mono_allelic_fiter_expr: hl.expr.BooleanExpression,
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
    if determine_cutoff_from_bin:
        snp_rf_cutoff, indel_rf_cutoff = aggregated_bin_ht.aggregate(
            [
                hl.agg.filter(
                    snv
                    & (aggregated_bin_ht.bin_id == bin_id)
                    & (aggregated_bin_ht.bin == cutoff),
                    hl.agg.min(aggregated_bin_ht.min_score),
                )
                for snv, cutoff in [
                    (aggregated_bin_ht.snv, snp_cutoff),
                    (~aggregated_bin_ht.snv, indel_cutoff),
                ]
            ]
        )
        snp_cutoff_global = hl.struct(bin=snp_cutoff, min_score=snp_rf_cutoff)
        indel_cutoff_global = hl.struct(
            bin=indel_cutoff, min_score=indel_rf_cutoff)

        logger.info(
            f"Using a SNP RF probability cutoff of {snp_rf_cutoff} and an indel RF probability cutoff of {indel_rf_cutoff}."
        )
    else:
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

    filters["InbreedingCoeff"] = hl.or_else(
        ht.InbreedingCoeff < inbreeding_coeff_cutoff, False
    )
    filters["AC0"] = ac0_filter_expr
    filters[
        "MonoAllelic"
    ] = mono_allelic_fiter_expr  # TODO: Do others agree that we should add this to gnomAD like we did for UKBB?

    # Fix annotations for release
    annotations_expr = {
        "rf_positive_label": hl.or_else(ht.tp, False),
        "rf_negative_label": ht.fail_hard_filters,
        "transmitted_singleton": hl.or_missing(
            ts_ac_filter_expr, ht.transmitted_singleton
        ),
        "rf_probability": ht.rf_probability["TP"],
    }
    if "feature_imputed" in ht.row:
        annotations_expr.update(
            {
                x: hl.or_missing(~ht.feature_imputed[x], ht[x])
                for x in [f for f in ht.row.feature_imputed]
            }
        )

    ht = ht.transmute(filters=add_filters_expr(
        filters=filters), **annotations_expr)

    ht = ht.annotate_globals(
        rf_snv_cutoff=snp_cutoff_global, rf_indel_cutoff=indel_cutoff_global
    )

    return ht


'''
def train_rf_model_local(
    ht: hl.Table,
    rf_features: List[str],
    tp_expr: hl.expr.BooleanExpression,
    fp_expr: hl.expr.BooleanExpression,
    fp_to_tp: float = 1.0,
    num_trees: int = 500,
    max_depth: int = 5,
    test_expr: hl.expr.BooleanExpression = False,
) -> Tuple[hl.Table, pyspark.ml.PipelineModel]:
    """
    Perform random forest (RF) training using a Table annotated with features and training data.
    .. note::
        This function uses `train_rf` and extends it by:
            - Adding an option to apply the resulting model to test variants which are withheld from training.
            - Uses a false positive (FP) to true positive (TP) ratio to determine what variants to use for RF training.
    The returned Table includes the following annotations:
        - rf_train: indicates if the variant was used for training of the RF model.
        - rf_label: indicates if the variant is a TP or FP.
        - rf_test: indicates if the variant was used in testing of the RF model.
        - features: global annotation of the features used for the RF model.
        - features_importance: global annotation of the importance of each feature in the model.
        - test_results: results from testing the model on variants defined by `test_expr`.
    :param ht: Table annotated with features for the RF model and the positive and negative training data.
    :param rf_features: List of column names to use as features in the RF training.
    :param tp_expr: TP training expression.
    :param fp_expr: FP training expression.
    :param fp_to_tp: Ratio of FPs to TPs for creating the RF model. If set to 0, all training examples are used.
    :param num_trees: Number of trees in the RF model.
    :param max_depth: Maxmimum tree depth in the RF model.
    :param test_expr: An expression specifying variants to hold out for testing and use for evaluation only.
    :return: Table with TP and FP training sets used in the RF training and the resulting RF model.
    """

    ht = ht.annotate(_tp=tp_expr, _fp=fp_expr, rf_test=test_expr)

    rf_ht = sample_training_examples(
        ht, tp_expr=ht._tp, fp_expr=ht._fp, fp_to_tp=fp_to_tp, test_expr=ht.rf_test
    )
    ht = ht.annotate(rf_train=rf_ht[ht.key].train,
                     rf_label=rf_ht[ht.key].label)

    summary = ht.group_by("_tp", "_fp", "rf_train", "rf_label", "rf_test").aggregate(
        n=hl.agg.count()
    )
    logger.info("Summary of TP/FP and RF training data:")
    summary.show(n=20)

    logger.info(
        "Training RF model:\nfeatures: {}\nnum_tree: {}\nmax_depth:{}".format(
            ",".join(rf_features), num_trees, max_depth
        )
    )

    rf_model = train_rf(
        ht.filter(ht.rf_train),
        features=rf_features,
        label="rf_label",
        num_trees=num_trees,
        max_depth=max_depth,
    )

    test_results = None
    if test_expr is not None:
        logger.info(f"Testing model on specified variants or intervals...")
        test_ht = ht.filter(hl.is_defined(ht.rf_label) & ht.rf_test)
        test_results = test_model(
            test_ht, rf_model, features=rf_features, label="rf_label"
        )

    print(test_results)

    features_importance = get_features_importance(rf_model)
    print("Test results:")
    print(test_results)
    print("Feature importance:")
    print(features_importance)
    print("RF featutes:")
    print(rf_features)
    ht = ht.select_globals(
        features_importance=features_importance,
        features=rf_features,
        # test_results=test_results
    )

    return ht.select("rf_train", "rf_label", "rf_test"), rf_model
'''

######################################
# main
########################################


def main(args):

    print("main table")
    ht = hl.read_table(
        f'{temp_dir}/ddd-elgh-ukbb/variant_qc/Sanger_table_for_RF_all_cols.ht')

    if args.train_rf:
        # ht = hl.read_table(
        #    f'{temp_dir}/ddd-elgh-ukbb/variant_qc/Sanger_table_for_RF_by_variant_type.ht')
        run_hash = str(uuid.uuid4())[:8]
        rf_runs = get_rf_runs(f'{tmp_dir}/ddd-elgh-ukbb/')
        while run_hash in rf_runs:
            run_hash = str(uuid.uuid4())[:8]
        ht_result, rf_model = train_rf(ht, args)
        print("Writing out ht_training data")
        ht_result = ht_result.checkpoint(
            get_rf(data="training", run_hash=run_hash).path, overwrite=True)
        # f'{tmp_dir}/ddd-elgh-ukbb/Sanger_RF_training_data.ht', overwrite=True)
        rf_runs[run_hash] = get_run_data(
            vqsr_training=False,
            transmitted_singletons=True,
            test_intervals=args.test_intervals,
            adj=True,
            features_importance=hl.eval(ht_result.features_importance),
            test_results=hl.eval(ht_result.test_results),
        )

        with hl.hadoop_open(f'{plot_dir}/ddd-elgh-ukbb/variant_qc/rf_runs.json', "w") as f:
            json.dump(rf_runs, f)
        pretty_print_runs(rf_runs)
        logger.info("Saving RF model")
        save_model(
            rf_model, get_rf(data="model", run_hash=run_hash), overwrite=True)
        # f'{tmp_dir}/ddd-elgh-ukbb/rf_model.model')
    else:
        run_hash = args.run_hash

    if args.apply_rf:

        logger.info(f"Applying RF model {run_hash}...")
        rf_model = load_model(get_rf(data="model", run_hash=run_hash))
        run_hash = args.run_hash
        # rf_model = hl.read_table(
        #    f'{temp_dir}/ddd-elgh-ukbb/variant_qc/models/{run_hash}/model.model')
        # ht = hl.read_table(
        #    f'{temp_dir}/ddd-elgh-ukbb/variant_qc/Sanger_cohorts_chr1-20-XY_sampleQC_FILTERED_FREQ_adj_inb.ht')
        ht = hl.read_table(
            f'{temp_dir}/ddd-elgh-ukbb/variant_qc/models/{run_hash}/training.ht')
        #ht = ht.annotate(rf_label=rf_model[ht.key].rf_label)
        ht = get_rf(data="training", run_hash=run_hash).ht()
        features = hl.eval(ht.features)
        ht = apply_rf_model(ht, rf_model, features, label=LABEL_COL)
        logger.info("Finished applying RF model")
        ht = ht.annotate_globals(rf_hash=run_hash)
        ht = ht.checkpoint(
            get_rf("rf_result_sanger_cohorts_new",
                   run_hash=run_hash).path, overwrite=True,
        )

        ht_summary = ht.group_by(
            "tp", "fp", TRAIN_COL, LABEL_COL, PREDICTION_COL
        ).aggregate(n=hl.agg.count())
        ht_summary.show(n=20)

    if args.finalize:
        run_hash = args.run_hash
        ht = hl.read_table(
            f'{temp_dir}/ddd-elgh-ukbb/variant_qc/models/{run_hash}/rf_result_ac_added.ht')
        # ht = create_grouped_bin_ht(
        #    model_id=run_hash, overwrite=True)
        freq_ht = hl.read_table(
            f'{temp_dir}/ddd-elgh-ukbb/variant_qc/Sanger_cohorts_chr1-20-XY_sampleQC_FILTERED_FREQ_adj.ht')
        freq = freq_ht[ht.key]

        # if not file_exists(
        #    get_score_quantile_bins(args.run_hash, aggregated=True).path
        # ):
        #    sys.exit(
        #        f"Could not find binned HT for RF  run {args.run_hash} (). Please run create_ranked_scores.py for that hash."
        #    )
        # aggregated_bin_ht = get_score_quantile_bins(ht, aggregated=True)
        print("created bin ht")

        ht = generate_final_rf_ht(
            ht,
            ac0_filter_expr=freq.freq[0].AC == 0,
            ts_ac_filter_expr=freq.freq[1].AC == 1,
            mono_allelic_fiter_expr=(freq.freq[1].AF == 1) | (
                freq.freq[1].AF == 0),
            snp_cutoff=args.snp_cutoff,
            indel_cutoff=args.indel_cutoff,
            determine_cutoff_from_bin=False,
            aggregated_bin_ht=bin_ht,
            bin_id=bin_ht.bin,
            inbreeding_coeff_cutoff=INBREEDING_COEFF_HARD_CUTOFF,
        )
        # This column is added by the RF module based on a 0.5 threshold which doesn't correspond to what we use
        # ht = ht.drop(ht[PREDICTION_COL])
        ht.write(f'{tmp_dir}/rf_final.ht', overwrite=True)


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
        "--list_rf_runs",
        help="Lists all previous RF runs, along with their hash, parameters and testing results.",
        action="store_true",
    )
    actions.add_argument(
        "--annotate_for_rf",
        help="Creates an annotated ht with features for RF",
        action="store_true",
    )
    actions.add_argument(
        "--train_rf", help="Trains RF model", action="store_true")
    actions.add_argument(
        "--apply_rf", help="Applies RF model to the data", action="store_true"
    )
    actions.add_argument(
        "--finalize", help="Write final RF model", action="store_true")

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
        default="chr20:1-1000000",
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
