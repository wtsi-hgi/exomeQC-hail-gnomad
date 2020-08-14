# sample qc pipeline based on gnomad qc pipeline
# for ddd-elgh-ukbb samples Sanger cohorts
# Pavlos Antoniou
# 14/08/2020
# HGI

import os
import hail as hl
import pandas as pd
import pyspark
import json
import sys
import re
from pathlib import Path
import logging
from typing import Any, Counter, List, Optional, Tuple, Union, Dict, Iterable

from bokeh.plotting import output_file, save, show
from sample_qc_v3.gnomad_methods.gnomad_filtering import compute_stratified_metrics_filter
from sample_qc_v3.gnomad_methods.gnomad_pipeline import get_qc_mt
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
    temp_dir = os.path.join(os.environ["HAIL_HOME"], "tmp")
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])

    cohorts_pop = hl.import_table(
        "s3a://DDD-ELGH-UKBB-exomes/ancestry/sanger_cohort_known_populations_ukbb_elgh_labels_updated.tsv", delimiter="\t").key_by('s')

    # drop cohorts
    # annotate with cohorts and populations from s3 table.
    # save matrixtable
    # mt = hl.read_matrix_table(
    #    f"{temp_dir}/ddd-elgh-ukbb/relatedness_ancestry/ddd-elgh-ukbb/chr1_chr20_XY_ldpruned.mt")
    #mt = mt.annotate_cols(original_pop=cohorts_pop[mt.s].known_population)
    ##mt = mt.annotate_cols(cohort=cohorts_pop[mt.s].cohort)
    #mt = mt.annotate_cols(known_pop=cohorts_pop[mt.s].known_population_updated)
    #mt = mt.annotate_cols(superpopulation=cohorts_pop[mt.s].superpopulation)
    #mt = mt.annotate_cols(gVCF=cohorts_pop[mt.s].gVCF_ID)
    # done the above on pca_RF jupyter notebook
    # mt.write(
    #    f"{tmp_dir}/ddd-elgh-ukbb/Sanger_chr1-20-XY_new_cohorts_split_multi_pops.mt", overwrite=True)

    # sex annotation gnomad v3 not available for dense matrixtable. stick with sex annotation from v2 already in matrixttable.

    mt = hl.read_matrix_table(
        f"{tmp_dir}/ddd-elgh-ukbb/Sanger_chr1-20-XY_new_cohorts_split_multi_pops.mt")
    # sample qc gnomad v3
    mt_qc_ready = get_qc_mt(mt)
    mt_qc_ready.write(
        f"{tmp_dir}/ddd-elgh-ukbb/Sanger_chr1-20-XY_new_cohorts_pops_QC_ready.mt")

    # ready for REF_ancestry AKT_overlap_updated_labels_August2020 and then pop_qc.py
