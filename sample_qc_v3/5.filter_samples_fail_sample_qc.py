# Pavlos Antoniou
# 16/09/2020
#  exclude samples that were outliers in population sample QC
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
from gnomad.resources.grch38 import gnomad
import argparse

tmp_dir = "hdfs://spark-master:9820/"
temp_dir = "file:///home/ubuntu/data/tmp"
plot_dir = "/home/ubuntu/data/tmp"

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


def main(args):
    mt = hl.read_matrix_table(
        args.matrixtable)
    samples_to_remove_filename = args.failed_samples_list

    samples_to_remove = hl.import_table(samples_to_remove_filename).key_by('s')
    mt_filtered = mt.filter_cols(hl.is_defined(
        samples_to_remove[mt.s]), keep=False)
    mt_filtered.write(
        f'{args.output_dir}/ddd-elgh-ukbb/Sanger_cohorts_chr1-7and20_split_sampleqc_filtered.mt')


if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])

    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory

    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])

    #####################################################################
    ###################### INPUT DATA  ##############################
    #####################################################################
    parser = argparse.ArgumentParser()
    # Read the matrixtable, chrX and chrY should be included
    input_params = parser.add_argument_group("Input parameters")
    input_params.add_argument(
        "--matrixtable",
        help="Full path of input matrixtable. Path format \"file:///home/ubuntu/data/tmp/path/to/.mt\"",
        default=f'{temp_dir}/ddd-elgh-ukbb/variant_qc/Sanger_cohorts_chr1-7and20_split.mt',
        type=str,
    )
    input_params.add_argument(
        "--failed_samples_list",
        help="Full path of tsv file with list of failed samples, one per line",
        default=f"{temp_dir}/ddd-elgh-ukbb/filtering/samples_failed_QC.tsv",
        type=str,
    )

    input_params.add_argument(
        "--output_dir",
        help="Full path of output folder to store results. Preferably hdfs or secure lustre",
        default=tmp_dir,
        type=str
    )

    args = parser.parse_args()
    main(args)
