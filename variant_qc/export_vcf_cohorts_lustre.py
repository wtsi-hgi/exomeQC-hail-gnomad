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

    mt = hl.read_matrix_table(
        f'{lustre_dir}/Sanger_cohorts_chr1-7and20_after_RF_final.mt')


    table_cohort = hl.import_table(
        f"{lustre_dir}/sanger_cohorts_corrected_ukbb_july_2020.tsv", delimiter="\t").key_by('s')

    mt = mt.annotate_cols(cohort=table_cohort[mt.s].cohort)
    df = pd.read_csv(
        f"{lustre_dir}/sanger_cohorts_corrected_ukbb_july_2020.tsv", sep="\t")
    cohorts_array = df.cohort.unique()

    mt = mt.annotate_rows(
        MAF_cohorts=hl.agg.group_by(mt.cohort,
                                    hl.min(hl.agg.call_stats(mt.GT, mt.alleles).AF))
    )
    mt = mt.annotate_rows(
        AN_cohorts=hl.agg.group_by(mt.cohort,
                                   hl.min(hl.agg.call_stats(mt.GT, mt.alleles).AN))
    )

    mt = mt.annotate_rows(
        AC_cohorts=hl.agg.group_by(mt.cohort,
                                   hl.min(hl.agg.call_stats(mt.GT, mt.alleles).AC))
    )

    mt = mt.annotate_rows(
        missingness_cohorts=hl.agg.group_by(mt.cohort, hl.min(
            (hl.agg.count_where(hl.is_missing(mt['GT']))) / mt.count_rows()*2))

    )

    mt = mt.annotate_rows(
        info=mt.info.annotate(cohort_names=mt.MAF_cohorts.keys())
    )
    mt = mt.annotate_rows(
        info=mt.info.annotate(MAF_cohorts_values=mt.MAF_cohorts.values())
    )

    mt = mt.annotate_rows(
        info=mt.info.annotate(AN_cohorts_values=mt.AN_cohorts.values())
    )

    mt = mt.annotate_rows(
        info=mt.info.annotate(AC_cohorts=mt.AC_cohorts.values())
    )

    mt = mt.annotate_rows(
        info=mt.info.annotate(
            missingness_cohorts_values=mt.missingness_cohorts.values())
    )

    mt = mt.checkpoint(
        f'{lustre_dir}/Sanger_WES_mt_with_stats.mt', overwrite=True)
    hl.export_vcf(
        mt, f'{lustre_dir}/Sanger_WES_chr1-7and20_after_RF_cohort_stats.vcf.bgz', parallel='separate_header')

    mt1 = mt.select_entries()
    mt_fin = mt1.filter_cols(mt1['s'] == 'sample')
    hl.export_vcf(
        mt_fin, f"{lustre_dir}/Sanger_WES_chr1-7and20_stats_for_VEP.vcf.bgz", parallel='separate_header')
