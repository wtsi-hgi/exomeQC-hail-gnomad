# Pavlos Antoniou
# 16/09/2020
#  trio matrixtable creation from fam file
import os
import hail as hl
import pandas as pd
import pyspark
import json
import sys
import re
from pathlib import Path
import logging
from typing import Any, Counter, List, Optional, Tuple, Union
from bokeh.plotting import output_file, save, show
from gnomad.resources.grch38 import gnomad
from gnomad.utils.annotations import unphase_call_expr, add_variant_type
from gnomad.variant_qc.random_forest import (
    apply_rf_model,
    load_model,
    median_impute_features,
    pretty_print_runs,
    save_model,
)

from hail import Table

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
    "AS_QD",
    "AS_MQRankSum",
    "AS_SOR",
    "AS_ReadPosRankSum",
]
TRUTH_DATA = ["hapmap", "omni", "mills", "kgp_phase1_hc"]
INBREEDING_COEFF_HARD_CUTOFF = -0.3


def split_info() -> hl.Table:
    """
    Generates an info table that splits multi-allelic sites from
    the multi-allelic info table.

    :return: Info table with split multi-allelics
    :rtype: Table
    """
    info_ht = get_info(split=False).ht()

    # Create split version
    info_ht = hl.split_multi(info_ht)

    # Index AS annotations
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(
            **{f: info_ht.info[f][info_ht.a_index - 1] for f in info_ht.info if f.startswith("AC") or (f.startswith("AS_") and not f == 'AS_SB_TABLE')},
            AS_SB_TABLE=info_ht.info.AS_SB_TABLE[0].extend(
                info_ht.info.AS_SB_TABLE[info_ht.a_index])
        ),
        lowqual=info_ht.lowqual[info_ht.a_index - 1]
    )
    return info_ht


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

    # ANNOTATION TABLES:
    truth_data_ht = hl.read_table(
        f'{temp_dir}/ddd-elgh-ukbb/variant_qc/truthset_table.ht')
    trio_stats_table = hl.read_table(
        f'{temp_dir}/ddd-elgh-ukbb/variant_qc/Sanger_cohorts_trios_stats.ht')
    #inbreeding_ht = hl.read_table(f'{temp_dir}/ddd-elgh-ukbb/variant_qc/Sanger_cohorts_inbreeding.ht')
    allele_data_ht = hl.read_table(
        f'{temp_dir}/ddd-elgh-ukbb/variant_qc/Sanger_cohorts_allele_data.ht')
    allele_counts_ht = hl.read_table(
        f'{temp_dir}/ddd-elgh-ukbb/variant_qc/Sanger_cohorts_qc_ac.ht')
    allele_counts_ht = allele_counts_ht.select(
        *['ac_qc_samples_raw', 'ac_qc_samples_adj'])
    inbreeding_ht = hl.read_table(
        f'{temp_dir}/ddd-elgh-ukbb/variant_qc/Sanger_cohorts_inbreeding.ht')
    group = "raw"
    mt = hl.read_matrix_table(
        f'{temp_dir}/ddd-elgh-ukbb/Sanger_cohorts_chr1to3-20_split.mt')

    # mt = mt.select_entries(
    #    GT=hl.unphased_diploid_gt_index_call(mt.GT.n_alt_alleles()))
    # mt = mt.annotate_rows(InbreedingCoeff=hl.or_missing(
    #   ~hl.is_nan(mt.info.InbreedingCoeff), mt.info.InbreedingCoeff))
    ht = mt.rows()
    ht = ht.transmute(**ht.info)
    ht = ht.select("FS", "MQ", "QD", "InbreedingCoeff", *INFO_FEATURES)

    trio_stats_ht = trio_stats_table.select(
        f"n_transmitted_{group}", f"ac_children_{group}"
    )

    ht = ht.annotate(
        **inbreeding_ht[ht.key],
        **trio_stats_ht[ht.key],
        **truth_data_ht[ht.key],
        **allele_data_ht[ht.key].allele_data,
        **allele_counts_ht[ht.key],
    )
    # Filter to only variants found in high quality samples or controls with no LowQual filter
    # ht = ht.filter(
    #    (ht[f"ac_children_{group}"] > 0)
    # )  # TODO: change to AS_lowqual for v3.1 or leave as is to be more consistent with v3.0? I will need to add this annotation if so

    ht = ht.select(
        "a_index",
        "was_split",
        *FEATURES,
        *TRUTH_DATA,
        **{
            "transmitted_singleton": (ht[f"n_transmitted_{group}"] == 1),
            "fail_hard_filters": (ht.QD < 2) | (ht.FS > 60) | (ht.MQ < 30),
        },

    )

    ht = ht.repartition(n_partitions, shuffle=False)
    ht = ht.checkpoint(
        f'{tmp_dir}/ddd-elgh-ukbb/Sanger_table_for_RF.ht', overwrite=True)
    ht = median_impute_features(ht, {"variant_type": ht.variant_type})
    ht = ht.checkpoint(
        f'{tmp_dir}/ddd-elgh-ukbb/Sanger_table_for_RF_by_variant_type.ht', overwrite=True)
