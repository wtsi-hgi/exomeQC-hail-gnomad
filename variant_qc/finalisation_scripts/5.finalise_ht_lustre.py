# Pavlos Antoniou
# 16/09/2020
#  trio matrixtable creation from fam file
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
from gnomad.resources.grch38 import gnomad
from gnomad.utils.annotations import unphase_call_expr, add_variant_type
from gnomad.variant_qc.pipeline import create_binned_ht, score_bin_agg, train_rf_model
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
#plot_dir = "/home/ubuntu/data/tmp"
plot_dir="/lustre/scratch123/teams/hgi/mercury/megaWES-variantqc"
lustre_dir = "file:///lustre/scratch123/teams/hgi/mercury/megaWES-variantqc"



######################################
# main
########################################

# annotate with RF probability score 1-
    # 'rf_probability': dict<str, float64>
    #  'rf_train': bool
    # 'rf_label': str
    # 'rf_test': bool
    # 'rf_probability': dict<str, float64> :  - rf_probability['TP']
    # 'rf_prediction': str
    # 'consequence': str
    # 'inheritance': str
    # 'rank': int64
    # 'score': float64

    # info rf_probability['TP']
    # FILTER column to PASS for ≤80 InDels ≤90 SNVs
   # ht = generate_final_rf_ht(
   #     ht,
   #     snp_cutoff=args.snp_cutoff,
   #     indel_cutoff=args.indel_cutoff,
   #     determine_cutoff_from_bin=False,
   #     bin_id=ht.bin,
   #     inbreeding_coeff_cutoff=INBREEDING_COEFF_HARD_CUTOFF,
   # )
    # This column is added by the RF module based on a 0.5 threshold which doesn't correspond to what we use
    #ht = ht.drop(ht[PREDICTION_COL])
def main(args):

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
    n_partitions = 500
    parser = argparse.ArgumentParser()

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
