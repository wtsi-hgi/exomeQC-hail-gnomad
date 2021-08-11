import os
import hail as hl
import pyspark
import json
import sys
import re
from pathlib import Path
import numpy as np
import pandas as pd
import cufflinks as cf
import plotly
import plotly.offline as py
import plotly.graph_objs as go
import plotly.express as px

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





table_cohort = "s3a://DDD-ELGH-UKBB-exomes/samples_cohorts.tsv"
intersection_bed = "s3a://DDD-ELGH-UKBB-exomes/agilentv5_UKB_intersection.bed"
union_bed = "s3a://DDD-ELGH-UKBB-exomes/agilentv5_UKB_union.bed"
agilent = "s3a://DDD-ELGH-UKBB-exomes/S04380110_Covered_nometadata.bed"


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
    cf.go_offline()  # required to use plotly offline (no account required).
    py.init_notebook_mode()
    #####################################################################
    ###################### INPUT DATA  ##############################
    #####################################################################
    CHROMOSOME = "chr1"

    #Read in annotated matrixtable for chromosome 1 with cohorts. 
    #merge with chrY

    #
    mt = hl.read_matrix_table(
        f"{temp_dir}/ddd-elgh-ukbb/{CHROMOSOME}.mt")

    mt_annotated = annotate_samples_with_cohort_info(mt, table_cohort)
    # mt_annotated.write(
    #    f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}_annotated.mt", overwrite=True)

    mt_annotated = mt_annotated.key_rows_by('locus').distinct_by_row(
    ).key_rows_by('locus', 'alleles')
    mt_split = hl.split_multi_hts(
        mt_annotated, keep_star=False, left_aligned=False)

    mt = mt_split.annotate_rows(
        Variant_Type=hl.cond((hl.is_snp(mt_split.alleles[0], mt_split.alleles[1])), "SNP",
                             hl.cond(
            hl.is_insertion(
                mt_split.alleles[0], mt_split.alleles[1]),
            "INDEL",
            hl.cond(hl.is_deletion(mt_split.alleles[0],
                                   mt_split.alleles[1]), "INDEL",
                    "Other"))))

    mt = mt.checkpoint(
        f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}-split-multi_cohorts.mt",  overwrite=True)
    print("Finished splitting and writing mt. ")

    agilent_table = hl.import_bed(
        agilent, reference_genome='GRCh38')
    mt_agilent = mt.filter_rows(hl.is_defined(agilent_table[mt.locus]))

    mt_agilent = hl.sample_qc(mt_agilent, name='sample_QC_Hail')
    pandadf1 = mt_agilent.cols().flatten()
    print("Outputting table of sample qc")
    pandadf1.export(
        f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}_agilent_sampleQC.tsv.bgz", header=True)

    mt = mt.checkpoint(
        f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}-sampleqc-unfiltered_annotated.mt", overwrite=True)
