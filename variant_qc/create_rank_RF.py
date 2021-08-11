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
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from bokeh.plotting import output_file, save, show
from gnomad.resources.grch38 import gnomad
from gnomad.utils.annotations import unphase_call_expr, add_variant_type
from gnomad.utils.annotations import annotate_adj, bi_allelic_expr

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


def create_rf_rank(data_type: str, run_hash: str) -> None:
    """
    Creates a ranked table for a RF run and writes it to its correct location in annotations.
    :param str data_type: One of 'exomes' or 'genomes'
    :param str run_hash: RF run hash
    :return: Nothing
    :rtype: None
    """
    logger.info(f"Creating rank file for {data_type} RF run {run_hash}")

    if not hl.hadoop_exists(f'gs://gnomad-tmp/gnomad_{data_type}_rf_{run_hash}.ht/_SUCCESS'):
        gnomad_ht = get_gnomad_annotations(data_type)
        ht = hl.read_table(rf_path(data_type, 'rf_result', run_hash=run_hash))
        ht = ht.annotate(**gnomad_ht[ht.key],
                         score=ht.rf_probability['TP'])

        # Write to temp location as result will be overwritten
        ht.write(f'gs://gnomad-tmp/gnomad_{data_type}_rf_{run_hash}.ht', overwrite=True)
    ht = hl.read_table(f'gs://gnomad-tmp/gnomad_{data_type}_rf_{run_hash}.ht')

    ht = add_rank(ht,
                  score_expr=1-ht.score,
                  subrank_expr={
                      'singleton_rank': ht.singleton,
                      'biallelic_rank': ~ht.was_split,
                      'biallelic_singleton_rank': ~ht.was_split & ht.singleton,
                      'adj_rank': ht.ac > 0,
                      'adj_biallelic_rank': ~ht.was_split & (ht.ac > 0),
                      'adj_singleton_rank': ht.singleton & (ht.ac > 0),
                      'adj_biallelic_singleton_rank': ~ht.was_split & ht.singleton & (ht.ac > 0)
                  }
                  )
    ht.write(rf_path(data_type, 'rf_result', run_hash=run_hash), overwrite=True)


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
    omni = f'{temp_dir}/ddd-elgh-ukbb/training_sets/1000G_omni2.5.hg38.ht'
    omni_ht = hl.read_table(omni)
    mills = f'{temp_dir}/ddd-elgh-ukbb/training_sets/Mills_and_1000G_gold_standard.indels.hg38.ht'
    mills_ht = hl.read_table(mills)
    thousand_genomes = f'{temp_dir}/ddd-elgh-ukbb/training_sets/1000G_phase1.snps.high_confidence.hg38.ht'
    thousand_genomes_ht = hl.read_table(thousand_genomes)
    hapmap = f'{temp_dir}/ddd-elgh-ukbb/training_sets/hapmap_3.3.hg38.ht'
    hapmap_ht = hl.read_table(hapmap)
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

    mt = hl.read_matrix_table(
        f'{temp_dir}/ddd-elgh-ukbb/filtering/Sanger_cohorts_chr1-20-XY_sampleQC_FILTERED.mt')
    mt = annotate_adj(mt)
    mt_freq = annotate_freq(mt)
    print("repartitioning:")
    #mt_freq = mt_freq.repartition(1000, shuffle=False)
    mt_freq = mt_freq.checkpoint(
        f'{tmp_dir}/Sanger_cohorts_chr1-20-XY_sampleQC_FILTERED_FREQ_adj.mt', overwrite=True)
    ht_freq = mt_freq.rows()
    ht_freq.describe()
    ht_freq.write(
        f'{tmp_dir}/Sanger_cohorts_chr1-20-XY_sampleQC_FILTERED_FREQ_adj.ht', overwrite=True)
