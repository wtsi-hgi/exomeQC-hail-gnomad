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

    # read matrixtable = remove the
    # mt = hl.read_matrix_table(
    #    f'{temp_dir}/ddd-elgh-ukbb/Sanger_cohorts_chr1-20-XY_new_cohorts_split_multi.mt')
    mt = hl.read_matrix_table(
        f'{temp_dir}/ddd-elgh-ukbb/variant_qc/Sanger_cohorts_chr1-7and20_split.mt')
    samples_to_remove_filename = f"{temp_dir}/ddd-elgh-ukbb/filtering/samples_failed_QC.tsv"

    samples_to_remove = hl.import_table(samples_to_remove_filename).key_by('s')
    mt_filtered = mt.filter_cols(hl.is_defined(
        samples_to_remove[mt.s]), keep=False)
    mt_filtered.write(
        f'{tmp_dir}/ddd-elgh-ukbb/Sanger_cohorts_chr1-7and20_split_sampleqc_filtered.mt')
