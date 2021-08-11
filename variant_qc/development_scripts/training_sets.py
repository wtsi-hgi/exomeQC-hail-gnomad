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
    mt2 = hl.read_matrix_table(
        f'{temp_dir}/ddd-elgh-ukbb/Sanger_cohorts_chr1-20-XY_new_cohorts_split_multi.mt')
    # mt2_hard_filter_fail = mt2.filter_rows(
    #    (mt2.info.QD < 2) | (mt2.info.FS > 60) | (mt2.info.MQ < 30))
    # mt2_hard_filter_fail = mt2_hard_filter_fail.checkpoint(
    #    f'{tmp_dir}/ddd-elgh-ukbb/Sanger_cohorts_variants_fail_hard_filter.mt', overwrite=True)
    # print(mt2.count())
    # print(mt2_hard_filter_fail.count())

    omni = f'{temp_dir}/ddd-elgh-ukbb/training_sets/1000G_omni2.5.hg38.ht'
    omni_ht = hl.read_table(omni)
    mills = f'{temp_dir}/ddd-elgh-ukbb/training_sets/Mills_and_1000G_gold_standard.indels.hg38.ht'
    mills_ht = hl.read_table(mills)
    thousand_genomes = f'{temp_dir}/ddd-elgh-ukbb/training_sets/1000G_phase1.snps.high_confidence.hg38.ht'
    thousand_genomes_ht = hl.read_table(thousand_genomes)

    mt_omni = mt2.filter_rows(hl.is_defined(omni_ht[mt2.row_key]), keep=True)
    mt_omni = mt_omni.checkpoint(
        f'{tmp_dir}/ddd-elgh-ukbb/Sanger_omni_TP.mt', overwrite=True)
    print(mt_omni.count())
    mt_mills = mt2.filter_rows(hl.is_defined(mills_ht[mt2.row_key]), keep=True)
    mt_mills = mt_mills.checkpoint(
        f'{tmp_dir}/ddd-elgh-ukbb/Sanger_mills_TP.mt', overwrite=True)
    print(mt_mills.count())
    mt_1000g = mt2.filter_rows(hl.is_defined(
        thousand_genomes_ht[mt2.row_key]), keep=True)
    mt_1000g = mt_1000g.checkpoint(
        f'{tmp_dir}/ddd-elgh-ukbb/Sanger_1000g_TP.mt', overwrite=True)
    print(mt_1000g.count())
