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

tmp_dir = "hdfs://spark-master:9820/"
temp_dir = "file:///home/ubuntu/data/tmp"
lustre_dir = "file:///lustre/scratch123/teams/hgi/mercury/pavlos-test"
plot_dir = "/home/ubuntu/data/tmp"

if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory

    hl.init(sc=sc, tmp_dir=lustre_dir, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])
    n_partitions = 500

    # 1. Import required files
    print("1. import required tables in hail for the project.")
    VQSLOD_snps = hl.import_table(f'{lustre_dir}/intervalwgs-qc/VQSLOD_snps.bgz',
                                  types={"Locus": "locus<GRCh38>", "VQSLOD": hl.tfloat64})
    VQSLOD_indels = hl.import_table(f'{lustre_dir}/intervalwgs-qc/VQSLOD_indels.bgz',
                                    types={"Locus": "locus<GRCh38>", "VQSLOD": hl.tfloat64})
    sample_QC_nonHail = hl.import_table(
        f'{lustre_dir}/intervalwgs-qc/INTERVAL_WGS_Sample_QC_04-09-2019.txt', impute=True)

    centromere_table = hl.import_bed(
        f'{lustre_dir}/intervalwgs-qc/Centromere_region_UCSC_GRCh38.bed', reference_genome='GRCh38', min_partitions=250)

    #####################################################################
    ###################### INPUT DATA  ##############################
    #####################################################################
    # Give chromosome as input to program with chr prefix i.e chr1, chr2, chr3 etc
    CHROMOSOME = "chr1"
    print(f"Reading {CHROMOSOME} mt")

    mt = hl.read_matrix_table(
        f'{lustre_dir}/chr1.mt')
    #mt = hl.read_matrix_table(f"{temp_dir}/matrixtables/{CHROMOSOME}.mt")

    print("Splitting mt and writing out split mt")
    mt_split = hl.split_multi_hts(mt, keep_star=False, left_aligned=False)

    mt_split = mt_split.checkpoint(
        f"{lustre_dir}/checkpoints/{CHROMOSOME}-split-multi_checkpoint.mt",  overwrite=True)
    print("Finished splitting and writing mt. ")

    #####################################################################
    ###################### UNFILTERED SAMPLE AND VARIANT QC #############
    #####################################################################

    print('Annotating rows with snp and indel info')
    mt = mt_split.annotate_rows(
        Variant_Type=hl.cond((hl.is_snp(mt_split.alleles[0], mt_split.alleles[1])), "SNP",
                             hl.cond(
            hl.is_insertion(
                mt_split.alleles[0], mt_split.alleles[1]),
            "INDEL",
            hl.cond(hl.is_deletion(mt_split.alleles[0],
                                   mt_split.alleles[1]), "INDEL",
                    "Other"))))

    # Unfiltered data summary stats:
    print("Finished annotating rows, annotating columns now")
    mt_sqc1_unfiltered = mt.annotate_cols(
        sample_QC_nonHail_unfiltered=sample_QC_nonHail.key_by("ID")[mt.s])
    mt_sqc2_unfiltered = hl.sample_qc(
        mt_sqc1_unfiltered, name='sample_QC_Hail')

    panda_df_unfiltered_table = mt_sqc2_unfiltered.cols().flatten()

    print("Outputting table of sample qc")
    panda_df_unfiltered_table.export(
        f"{lustre_dir}/output-tables/{CHROMOSOME}_sampleQC_unfiltered.tsv.bgz", header=True)

    # Variant QC
    mt2 = hl.variant_qc(mt_sqc2_unfiltered, name='variant_QC_Hail')

    print('Exporting variant qc pandas table to disk')
    mt_rows = mt2.rows()
    mt_rows.select(mt_rows.variant_QC_Hail).flatten().export(f"{lustre_dir}/output-tables/{CHROMOSOME}_variantQC_unfiltered.tsv.bgz",
                                                             header=True)

    #####################################################################
    ###################### START FILTERING ##############################
    #####################################################################
    mt = mt2
    # VQSR filtering
    print("Annotation VQSLOD snp and indel scores")
    mt = mt.annotate_rows(VQSLOD_SNP=VQSLOD_snps.key_by("Locus")[mt.locus])
    mt = mt.annotate_rows(VQSLOD_INDEL=VQSLOD_indels.key_by("Locus")[mt.locus])

    print("Filtering on VQSLOD scores")
    vqslod_filtered = (
        hl.case()
        .when((mt.Variant_Type == "SNP"), (mt.VQSLOD_SNP.VQSLOD >= thresholds['intervalwgs']['snp_vqsr_threshold']))
        .when((mt.Variant_Type == "INDEL"), (mt.VQSLOD_INDEL.VQSLOD >= thresholds['intervalwgs']['indel_vqsr_threshold']))
        .default(False)  # remove everything else
    )

    mt_vqslod_filtered = mt.filter_rows(vqslod_filtered)
    print("Finished filtering. Writing out matrixtable...")

    mt_vqslod_filtered = mt_vqslod_filtered.checkpoint(
        f"{lustre_dir}/checkpoints/{CHROMOSOME}_vqslod_filtered_checkpoint.mt", overwrite=True)
    print("Finished writing vqslod filtered matrixtable")

    # CENTROMERE FILTERING
    print("Centromere filtering")
    mt_vqslod_filtered_WO_centromere = mt_vqslod_filtered.filter_rows(
        hl.is_defined(centromere_table[mt_vqslod_filtered.locus]), keep=False)

    mt_sqc1 = mt_vqslod_filtered_WO_centromere.annotate_cols(
        sample_QC_nonHail=sample_QC_nonHail.key_by("ID")[mt_vqslod_filtered_WO_centromere.s])

    # Sample QC filtering
    print("Filtering on sample qc")
    mt_sqc1_filtered = mt_sqc1.filter_cols(
        (mt_sqc1.sample_QC_nonHail["PASS.Depth"] == 1) &
        (mt_sqc1.sample_QC_nonHail["PASS.ID"] == 1) &
        (mt_sqc1.sample_QC_nonHail["PASS.Median.FreeMix"] == 1) &
        (mt_sqc1.sample_QC_nonHail["PASS.NRD"] == 1) &
        (mt_sqc1.sample_QC_nonHail["PASS.SampleSwap"] == 1) &
        (mt_sqc1.sample_QC_nonHail["PASS.Sex"] == 1) &
        (mt_sqc1.sample_QC_nonHail["PASS.DUP"] == 1)
    )
    print("Writing out filtered sample qc checkpoint")
    mt_sqc2 = hl.sample_qc(mt_sqc1_filtered, name='sample_QC_Hail')

    mt_sqc2 = mt_sqc2.checkpoint(f"{lustre_dir}/checkpoints/{CHROMOSOME}_sampleQC_filtered_checkpoint.mt",
                                 overwrite=True)

    filter_sampleqc_table = mt_sqc2.cols().flatten()
    filter_sampleqc_table.export(
        f"{lustre_dir}/output-tables/{CHROMOSOME}-sampleQC_filtered.tsv.bgz", header=True)

    # ADD allele bias script here

    initial_geno = mt_sqc2.aggregate_entries(
        hl.agg.fraction(hl.is_defined(mt_sqc2.GT)))

    print(f'Defined genotypes: {initial_geno * 100:.2f}%.')
    initial_missing = 100 - (initial_geno * 100)
    print(f'Initial missing genotype: {initial_missing:.2f}%.')

    # Remove initial missing genotypes
    mt_sqc2 = mt_sqc2.filter_entries(hl.is_defined(mt_sqc2.GT))

    ab = mt_sqc2.AD[1] / hl.sum(mt_sqc2.AD)

    filter_condition_ab = (
        hl.case(missing_false=True)
        .when(mt_sqc2.GT.is_hom_ref(), ab > 0.1)
        .when(mt_sqc2.GT.is_het(), (ab < 0.20) | (ab > 0.80))
        .when(mt_sqc2.GT.is_hom_var(), ab < 0.9)
        .default(False)  # remove everything else
    )

    fraction_filtered = mt_sqc2.aggregate_entries(
        hl.agg.fraction(filter_condition_ab))
    print(
        f'Filtering {fraction_filtered * 100:.2f}% entries out of downstream analysis.')

    print(
        f'Total Filtering {(fraction_filtered * 100) + initial_missing:.2f}% entries including initial missing data out of downstream analysis.')

    Total_geno = mt_sqc2.aggregate_entries(
        hl.agg.count_where(hl.is_defined(mt_sqc2.GT)))
    print('Total genotypes: ' + str(Total_geno))

    mt_sqc2_GT = mt_sqc2.filter_entries(filter_condition_ab, keep=False)

    ############################

    pro_AD_DP = hl.sum(mt_sqc2_GT.AD) / mt_sqc2_GT.DP

    Other_filters = (
        hl.case(missing_false=True)
        .when(hl.is_defined(mt_sqc2_GT.GT), (mt_sqc2_GT.DP > 100) | (pro_AD_DP < 0.9))
        .default(False)  # remove everything else
    )

    mt_sqc3 = mt_sqc2_GT.filter_entries(Other_filters, keep=False)

    # Saving on s3
    #mt_sqc2 = mt_sqc2_GT.checkpoint('s3a://interval-wgs/checkpoints_new/chr20_sampleQC_step1_filtered_allele_balance_checkpoint.mt', _read_if_exists = True)

    # Saving on local volume
    mt_sqc3 = mt_sqc3.checkpoint(f"{lustre_dir}/checkpoints/{CHROMOSOME}_sampleQC_step1_filtered_GenotypeQC_checkpoint.mt",
                                 overwrite=True)

    ######################
    # VARIANT qc
    ######################
    print("Variant qc:")
    mt_sqc_vqc = hl.variant_qc(mt_sqc3, name='variant_QC_Hail')
    mt_sqc_vqc_filtered = mt_sqc_vqc.filter_rows(
        (mt_sqc_vqc.variant_QC_Hail.call_rate >= 0.95) &
        (mt_sqc_vqc.variant_QC_Hail.p_value_hwe >= 10 ** -6) &
        (mt_sqc_vqc.variant_QC_Hail.gq_stats.mean >= 20) &
        (mt_sqc_vqc.variant_QC_Hail.AC[1] >= 1)
    )

    #####################################################################
    ###################### FINAL QC AFTER FILTERING  ####################
    #####################################################################

    fields_to_drop = ['variant_QC_Hail', 'sample_QC_Hail']

    mt1 = mt_sqc_vqc_filtered.drop(*fields_to_drop)

    mt2 = hl.sample_qc(mt1, name='sample_QC_Hail')
    mt3 = hl.variant_qc(mt2, name='variant_QC_Hail')

    mt3 = mt3.checkpoint(
        f"{lustre_dir}/checkpoints/{CHROMOSOME}-full-sampleqc-variantqc-filtered-FINAL.mt", overwrite=True)
    mt3_cols = mt3.cols()
    mt3_cols.flatten().export(
        f"{lustre_dir}/output-tables/{CHROMOSOME}-sampleQC_filtered_FINAL.tsv.bgz", header=True)

    mt3_rows = mt3.rows()
    mt3_rows.select(mt3_rows.variant_QC_Hail).flatten().export(
        f"{lustre_dir}/output-tables/{CHROMOSOME}-variantQC_filtered_FINAL.tsv.bgz", header=True)
