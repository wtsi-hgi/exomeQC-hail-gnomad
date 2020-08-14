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
from gnomad_ancestry import pc_project, run_pca_with_relateds, assign_population_pcs

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

    bed_to_exclude_pca = hl.import_bed(
        f"{temp_dir}/1000g/price_high_ld.bed.txt", reference_genome='GRCh38')

    cohorts_pop = hl.import_table(
        "s3a://DDD-ELGH-UKBB-exomes/ancestry/sanger_cohort_known_populations_ukbb_elgh_labels.tsv", delimiter="\t").key_by('s')

    pca_scores = hl.read_table(
        f"{temp_dir}/ddd-elgh-ukbb/elgh_labels/pop_assignments_test.ht")
    # pca_loadings = hl.read_table(f"{temp_dir}/ddd-elgh-ukbb/pca_loadings.ht")
    logger.info("assign population pcs")
   # population_assignment_table = assign_population_pcs(
    #    pca_scores, pca_loadings, known_col="known_pop")

    pop_ht, pop_clf = assign_population_pcs(
        pca_scores, pca_scores.pca_scores, known_col="known_pop", n_estimators=100, prop_train=0.8, min_prob=0.5)
    pop_ht.write(
        f"{tmp_dir}/ddd-elgh-ukbb/pop_assignments_test_minprob_0.5.ht", overwrite=True)
    pop_ht.export(
        f"{temp_dir}/ddd-elgh-ukbb/pop_assignments_elgh_labels.tsv.gz")
    #filename = f"{temp_dir}/ddd-elgh-ukbb/RF_model.pkl"
    #pickle.dump(pop_clf, open(filename, 'wb'))
