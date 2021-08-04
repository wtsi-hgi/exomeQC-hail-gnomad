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
import argparse
from typing import Any, Counter, List, Optional, Tuple, Union


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

def main():

    print("main")

    run_hash = "91ba5f38"
    ht=hl.read_table(f'{lustre_dir}/variant_qc/models/{run_hash}_score_binning.ht')

    mt = hl.read_matrix_table(
        f'{lustre_dir}/MegaWESSanger_cohorts_sampleQC_filtered.mt')
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
    mt = mt.annotate_rows(
        info=mt.info.annotate(bin=ht[mt.row_key].bin)
    )
    
    filter_column_annotation = (
        hl.case()
        .when(((mt.Variant_Type == "SNP") & (mt.info.bin <= 91)), "PASS")
        .when(((mt.Variant_Type == "INDEL") & (mt.info.bin <= 0.86)), "PASS")
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
    print(f'Failed:{mt_fail2.count()}')
    print(f'Passed:{mt_pass.count()}')

    mt2 = mt2.checkpoint(
        f'{lustre_dir}/variant_qc/megaWES_final_after_RF_{run_hash}.mt', overwrite=True)
   
    chroms=[*range(1,23),"X","Y"]
    chromosomes=["chr"+ str(chr) for chr in chroms]
    for chromosome in chromosomes:
        print(chromosome)
        mt=mt2.filter_rows(mt2.locus.contig==chromosome)
        mt.write(f'{lustre_dir}/final_matrixtables_VCFs/{chromosome}_after_RF_{run_hash}.mt',overwrite=True)
        hl.export_vcf(
        mt, f'{lustre_dir}/final_matrixtables_VCFs/VCFs/{chromosome}_after_RF_{run_hash}.vcf.bgz',parallel='separate_header')


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
    
   
    main()
