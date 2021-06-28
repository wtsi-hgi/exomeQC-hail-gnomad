# Pavlos Antoniou
# 28/06/2021
#  Generate files required for RF part 1
import os
import hail as hl
import pandas as pd
import pyspark
import json
import sys
import re
from pathlib import Path
import logging
import argparse
from typing import Any, Counter, List, Optional, Tuple, Union
from gnomad.resources.grch38 import gnomad
from gnomad.utils.annotations import annotate_adj
from gnomad.utils.annotations import unphase_call_expr, add_variant_type
from gnomad.utils.annotations import annotate_adj, bi_allelic_expr, bi_allelic_site_inbreeding_expr
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.sample_qc.relatedness import (
    SIBLINGS,
    generate_sib_stats_expr,
    generate_trio_stats_expr,
)
from hail import Table
from helper_functions import *

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

LABEL_COL = "rf_label"
TRAIN_COL = "rf_train"
PREDICTION_COL = "rf_prediction"

TRUTH_DATA = ["hapmap", "omni", "mills", "kgp_phase1_hc"]
INBREEDING_COEFF_HARD_CUTOFF = -0.3

################################
# Define the hail persistent storage directory
tmp_dir = "hdfs://spark-master:9820/"
temp_dir = "file:///home/ubuntu/data/tmp"
lustre_dir = "file:///lustre/scratch123/teams/hgi/mercury/megaWES-variantqc"
######################################



def main(args):
    group = "raw"
 
    mt = hl.read_matrix_table(
        args.matrixtable)

    # Truthset
    mt = hl.variant_qc(mt)

    truthset_ht = get_truth_ht( args.omni, args.mills, args.thousand_genomes, args.hapmap)
    truthset_ht.write(
        f'{args.output_dir}/variant_qc/truthset_table.ht', overwrite=True)
    # Trio data
    # trio annotation:
    logger.info("Trio annotation and writing trios_adj.mt")
    mt_adj = annotate_adj(mt)
    fam = args.trio_fam
    pedigree = hl.Pedigree.read(fam)
    trio_dataset = hl.trio_matrix(mt_adj, pedigree, complete_trios=True)
    trio_dataset.checkpoint(
        f'{args.output_dir}/variant_qc/MegaWES_trios_adj.mt', overwrite=True)
    logger.info("Trio stats and writing MegaWes_stats.ht")
    trio_stats_ht = generate_trio_stats(
        trio_dataset, autosomes_only=True, bi_allelic_only=True)
    trio_stats_ht.write(
        f'{args.output_dir}/variant_qc/MegaWES_stats.ht', overwrite=True)

    # inbreeding ht
    mt_inbreeding = mt.annotate_rows(
        InbreedingCoeff=bi_allelic_site_inbreeding_expr(mt.GT))

    mt = mt.key_rows_by('locus').distinct_by_row(
    ).key_rows_by('locus', 'alleles')

    ht_inbreeding = mt_inbreeding.rows()

    # allele data and qc_ac ht
    allele_data_ht = generate_allele_data(mt)

    qc_ac_ht = generate_ac(mt, fam)

    logger.info("Writing tables for inbreeding, allele counts")
    ht_inbreeding.write(
        f'{args.output_dir}/variant_qc/MegaWES_inbreeding_new.ht', overwrite=True)
    qc_ac_ht.write(
        f'{args.output_dir}/variant_qc/MegaWES_qc_ac_new.ht', overwrite=True)
    allele_data_ht.write(
        f'{args.output_dir}/variant_qc/MegaWES_allele_data_new.ht', overwrite=True)


    # Trio matrix table
    logger.info("Split multi allelic variants and write mt")
    mt = hl.split_multi_hts(
        mt, keep_star=False, left_aligned=False, permit_shuffle=True)
    mt=mt.checkpoint(f'{args.output_dir}/variant_qc/MegaWESSanger_cohorts_sampleQC_filtered_split.mt', overwrite=True)
    fam = args.trio_fam
    pedigree = hl.Pedigree.read(fam)
    logger.info("Trio matrixtable generation:")
    trio_dataset = hl.trio_matrix(mt, pedigree, complete_trios=True)
    trio_dataset.write(
        f'{args.output_dir}/variant_qc/MegaWES_trio_table.mt', overwrite=True)

    # Family stats
    logger.info("Family stats")
    (ht1, famstats_ht) = generate_family_stats(mt, fam)
    print("Writing mt and family stats_ht")
    ht1.write(f'{args.output_dir}/variant_qc/MegaWES_family_stats.ht',
              overwrite=True)

    mt = mt.annotate_rows(family_stats=ht1[mt.row_key].family_stats)
    mt=mt.checkpoint(f'{args.output_dir}/variant_qc/MegaWES_family_stats.mt', overwrite=True)

    #Family stats with Allele Frequencies from gnomad
    logger.info("Family stats with gnomad AF")
    priors = hl.read_table(args.priors)
    mt = mt.annotate_rows(gnomad_maf=priors[mt.row_key].maf)
    mt = mt.checkpoint(
        f'{lustre_dir}/variant_qc/MegaWES_family_stats_gnomad_AF.mt', overwrite=True)

    logger.info("De novo table cration")
    #De novo table
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
        default=f'{lustre_dir}/MegaWESSanger_cohorts_sampleQC_filtered.mt',
        type=str,
    )
    input_params.add_argument(
        "--omni",
        help="Full path of Omni truset hail table GRCh38",
        default=f'{lustre_dir}/training_sets/1000G_omni2.5.hg38.ht',
        type=str,
    )
    input_params.add_argument(
        "--mills",
        help="Full path of Mills_1000G_indels truset hail table GRCh38",
        default=f'{lustre_dir}/training_sets/Mills_and_1000G_gold_standard.indels.hg38.ht',
        type=str,
    )
    input_params.add_argument(
        "--thousand_genomes",
        help="Full path of 1000 genomes project truset hail table GRCh38",
        default=f'{lustre_dir}/training_sets/1000G_phase1.snps.high_confidence.hg38.ht',
        type=str,
    )
    input_params.add_argument(
        "--hapmap",
        help="Full path of Hapmap truset hail table GRCh38",
        default=f'{lustre_dir}/training_sets/hapmap_3.3.hg38.ht',
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
