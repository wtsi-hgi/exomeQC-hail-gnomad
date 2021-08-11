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
from typing import Any, Counter, List, Optional, Tuple, Union, Dict
from bokeh.plotting import output_file, save, show
from gnomad.resources.grch38 import gnomad
from gnomad.utils.annotations import annotate_adj
from gnomad.utils.annotations import unphase_call_expr, add_variant_type
from gnomad.utils.annotations import annotate_adj, bi_allelic_expr, bi_allelic_site_inbreeding_expr
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.filtering import filter_to_adj
from gnomad.sample_qc.relatedness import (
    SIBLINGS,
    generate_sib_stats_expr,
    generate_trio_stats_expr,
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
    "AS_pab_max",
    "AS_MQRankSum",
    "AS_SOR",
    "AS_ReadPosRankSum",
]
TRUTH_DATA = ["hapmap", "omni", "mills", "kgp_phase1_hc"]
INBREEDING_COEFF_HARD_CUTOFF = -0.3


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

    ################################

    #################################

    # mt_trios = hl.read_matrix_table(
    #    f'{temp_dir}/ddd-elgh-ukbb/variant_qc/mt_trios_adj.mt')
    run_hash = "91b132aa"
    ht = hl.read_table(
        f'{temp_dir}/ddd-elgh-ukbb/variant_qc/models/{run_hash}/{run_hash}_rf_result_sanger_cohorts_DENOVO_family_stats_SYNONYMOUS.ht')

    # mt_trios = mt_trios.annotate_rows(
    #    consequence=ht[mt_trios.row_key].consequence)

    # mt_trios = mt_trios.checkpoint(
    #    f'{tmp_dir}/sanger_cohorts_trios_consequence.mt', overwrite=True)
    # mt_filtered = mt_trios.filter_rows((mt_trios.info.AC[0] <= 2) & (
    #    mt_trios.consequence == "synonymous_variant"))
    mt_filtered = hl.read_matrix_table(
        f'{temp_dir}/ddd-elgh-ukbb/variant_qc/sanger_cohorts_AC_synonymous_filtered.mt')
    # mt_filtered = mt_filtered.checkpoint(
    #    f'{tmp_dir}/sanger_cohorts_AC_synonymous_filtered.mt', overwrite=True)

    mt_trans = mt_filtered.filter_entries(mt_filtered.info.AC[0] == 2)
    mt_untrans = mt_filtered.filter_entries(mt_filtered.info.AC[0] == 1)
    mt_trans = mt_trans.group_cols_by(mt_trans.id).aggregate(transmitted_singletons_count=hl.agg.count_where((mt_trans.info.AC[0] == 2) &
                                                                                                             (mt_trans.proband_entry.GT.is_het_ref()) &
                                                                                                             (mt_trans.father_entry.GT.is_het_ref()) |
                                                                                                             (mt_trans.mother_entry.GT.is_het_ref())))
    mt_untrans = (mt_untrans.group_cols_by(mt_untrans.id).aggregate(
        untransmitted_singletons_count=hl.agg.count_where((mt_untrans.info.AC[0] == 1) &
                                                          (mt_untrans.proband_entry.GT.is_hom_ref()) &
                                                          (mt_untrans.father_entry.GT.is_het_ref()) |
                                                          (mt_untrans.mother_entry.GT.is_het_ref()))))

    ht_entries = mt_trans.entries().key_by('locus', 'alleles')
    ht1 = ht.annotate(
        transmitted_singletons=ht_entries[ht.key].transmitted_singletons_count)
    ht1 = ht1.annotate(family_id=ht_entries[ht1.key].id)

    ht_untrans_entries = mt_untrans.entries().key_by('locus', 'alleles')
    ht2 = ht1.annotate(
        untransmitted_singletons=ht_untrans_entries[ht1.key].untransmitted_singletons_count)
    print(ht2.describe())
    print(ht2.transmitted_singletons.show())
    print(ht2.untransmitted_singletons.show())
    ht2.write(
        f'{tmp_dir}/ddd-elgh-ukbb/{run_hash}/{run_hash}_rf_result_transmitted_singletons.ht', overwrite=True)
