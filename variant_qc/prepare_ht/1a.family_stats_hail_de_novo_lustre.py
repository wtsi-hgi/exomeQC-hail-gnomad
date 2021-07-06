# Pavlos Antoniou
# 16/09/2020
#  trio matrixtable creation from fam file
from hail import Table
import os
import hail as hl
import pandas as pd
import pyspark
import json
import sys
import re
import argparse
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
################################
# Define the hail persistent storage directory
tmp_dir = "hdfs://spark-master:9820/"
temp_dir = "file:///home/ubuntu/data/tmp"
plot_dir = "/home/ubuntu/data/tmp"
lustre_dir = "file:///lustre/scratch123/teams/hgi/mercury/megaWES-variantqc"
######################################

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



def main(args):
    ################################

    truthset_table = hl.read_table(args.truthset_table)
    #################################
    group = "raw"

    mt = hl.read_matrix_table(args.matrixtable)
    mt = hl.split_multi_hts(
        mt, keep_star=False, left_aligned=False, permit_shuffle=True)
    mt=mt.checkpoint(f'{args.output_dir}/variant_qc/MegaWESSanger_cohorts_sampleQC_filtered_autosomes_split.mt', overwrite=True)
    fam = args.trio_fam
    pedigree = hl.Pedigree.read(fam)
    trio_dataset = hl.trio_matrix(mt, pedigree, complete_trios=True)
    trio_dataset.write(
        f'{args.output_dir}/MegaWES_trio_table.mt', overwrite=True)

    (ht1, famstats_ht) = generate_family_stats(mt, fam)
    print("Writing mt and family stats_ht")
    ht1.write(f'{args.output_dir}/MegaWES_family_stats.ht',
              overwrite=True)
    #famstats_ht.write(
    #    f'{args.output_dir}/MegaWES_family_stats.ht', overwrite=True)
    mt = mt.annotate_rows(family_stats=ht1[mt.row_key].family_stats)
    mt=mt.checkpoint(f'{args.output_dir}/MegaWES_family_stats.mt', overwrite=True)
    #(mt1, famstats_ht) = generate_family_stats(mt, fam)
    #print("Writing mt and family stats_ht")
    #mt1.write(f'{tmp_dir}/Sanger_cohorts_family_stats.mt', overwrite=True)
    # famstats_ht.write(
    #    f'{tmp_dir}/Sanger_cohorts_family_stats.ht', overwrite=True)
    #mt = mt.annotate_rows(family_stats=ht[mt.row_key].family_stats)
    #mt.write(f'{tmp_dir}/Sanger_cohorts_family_stats.mt', overwrite=True)

    priors = hl.read_table(args.priors)
    mt = mt.annotate_rows(gnomad_maf=priors[mt.row_key].maf)
    mt = mt.checkpoint(
        f'{lustre_dir}/variant_qc/MegaWES_family_stats_gnomad_AF.mt', overwrite=True)
    #mt = hl.split_multi_hts(mt, keep_star=False, left_aligned=False, permit_shuffle=True)
    
    de_novo_table = hl.de_novo(
        mt, pedigree, mt.gnomad_maf)

    de_novo_table = de_novo_table.key_by(
        'locus', 'alleles').collect_by_key('de_novo_data')
    de_novo_table.write(
        f'{args.output_dir}/variant_qc/MegaWES_denovo_table.ht', overwrite=True)


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

    #####################################################################
    ###################### INPUT DATA  ##############################
    #####################################################################
    parser = argparse.ArgumentParser()
    # Read the matrixtable
    input_params = parser.add_argument_group("Input parameters")
    input_params.add_argument(
        "--matrixtable",
        help="Full path of input matrixtable. Path format \"file:///home/ubuntu/data/tmp/path/to/.mt\"",
        default=f'{lustre_dir}/variant_qc/MegaWESSanger_cohorts_sampleQC_filtered_autosomes.mt',
        type=str,
    )
    input_params.add_argument(
        "--truthset_table",
        help="Full path of the truthset table created in variant qc step 1.",
        default=f'{lustre_dir}/variant_qc/truthset.ht',
        type=str,
    )
    input_params.add_argument(
        "--trio_fam",
        help="Full path of trio fam file for cohort",
        default=f"{lustre_dir}/trios/DDD_trios.fam",
        type=str,
    )
    input_params.add_argument(
        "--priors",
        help="Full path of prior AF for gnomad cohort",
        default=f'{lustre_dir}/gnomad_v3-0_AF.ht',
        type=str,
    )
    input_params.add_argument(
        "--output_dir",
        help="Full path of output folder to store results. Preferably hdfs or secure lustre",
        default=lustre_dir,
        type=str
    )

    args = parser.parse_args()
    main(args)
