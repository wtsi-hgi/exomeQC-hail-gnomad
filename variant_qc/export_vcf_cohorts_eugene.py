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
plot_dir = "/home/ubuntu/data/tmp"

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

    mt = hl.read_matrix_table(
        f'{temp_dir}/ddd-elgh-ukbb/variant_qc/Sanger_cohorts_chr1-7and20_after_RF_final.mt')

    table_cohort = hl.import_table(
        f"{temp_dir}/ddd-elgh-ukbb/sanger_cohorts_corrected_ukbb_july_2020.tsv", delimiter="\t").key_by('s')

    mt_result = mt.annotate_cols(cohort=table_cohort[mt.s].cohort)
    df = pd.read_csv(
        f"{temp_dir}/ddd-elgh-ukbb/sanger_cohorts_corrected_ukbb_july_2020.tsv", sep="\t")
    cohorts_array = df.cohort.unique()

    for cohort in cohorts_array:
        mt_cohort = mt_result.filter_cols(mt_result['cohort'] == cohort)
        numofsamples = mt_cohort.aggregate_cols(
            hl.agg.count_where(hl.is_defined('s')))
        print(numofsamples)
        n_called = hl.agg.count_where(hl.is_defined(mt_cohort.GT))
        missingness = (n_called / numofsamples) * 2
        mt_cohort = mt_cohort.annotate_rows(
            call_stats=hl.agg.call_stats(mt_cohort.GT, mt_cohort.alleles))
        mt_cohort = mt_cohort.annotate_rows(AC=mt_cohort.call_stats.AC)
        mt_cohort = mt_cohort.annotate_rows(AN=mt_cohort.call_stats.AN)
        mt_cohort = mt_cohort.annotate_rows(
            maf=hl.float64(mt_cohort.final_AC[1]/mt_cohort.final_AN))
        print(mt_cohort.AC)
        print(mt_cohort.AN)
        print(mt_cohort.maf)
