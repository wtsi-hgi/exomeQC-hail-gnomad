import os
import hail as hl
import pyspark
import json
import sys
import re
from pathlib import Path

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

partitions = 1000
tmp_dir = "hdfs://spark-master:9820/"
temp_dir = "file:///home/ubuntu/data/tmp"
lustre_dir = "file:///lustre/scratch123/teams/hgi/mercury/megaWES-variantqc"

plot_dir = "/home/ubuntu/data/tmp"

CHROMOSOMES = ["chr2",
               "chr3",
               "chr4",
               "chr5",
               "chr6",
               "chr7",
               "chr8",
               "chr9",
               "chr10",
               "chr11",
               "chr12",
               "chr13",
               "chr14",
               "chr15",
               "chr16",
               "chr17",
               "chr18",
               "chr19",
               "chr20",
               "chr21",
               "chr22",
               "chrX",
               "chrY"
               ]


def annotate_samples_with_cohort_info(mt: hl.MatrixTable, cohort_file) -> hl.MatrixTable:
    '''

    :param mt: matrixtable with cohort samples and variants
    :param cohort_file: a txt file with no header line and 2 columns, 1st: for sampleID;  2nd: cohortname; example:/lustre/scratch115/projects/autozyg/new_autozyg_DDD_callset.April2019/sample_list.after_QC.ELGH_BiB_Birm_controls_only.with_cohort_labels.txt
    :return: matrixtable with new column annotation
    '''
    # import the tab delimited file. Note that it is important for joins of tables to have defined keys in the hail tables
    table_cohort = hl.import_table(cohort_file, key='sample')
    # annotate the samples with a new attribute called cohort:
    mt_result = mt.annotate_cols(cohort=table_cohort[mt.s].cohort)
    return mt_result


table_cohort = "s3a://DDD-ELGH-UKBB-exomes/samples_cohorts.tsv"

if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory

    hl.init(sc=sc, tmp_dir=lustre_dir, local_tmpdir=lustre_dir, min_block_size=256, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])

    #####################################################################
    ###################### INPUT DATA  ##############################
    #####################################################################
    PARTITIONS=40000
    

    print("Reading all matrixtables and joining them")
    mt = hl.read_matrix_table(f"{lustre_dir}/chr1.mt")
    for chromosome in CHROMOSOMES:
        mt2 = hl.read_matrix_table(f"{lustre_dir}/{chromosome}.mt")
        mt = mt.union_rows(mt2)

    if mt.n_partitions() > PARTITIONS:
        print("Performing naive coalescing")
        mt = mt.naive_coalesce(PARTITIONS)
    print("Now annotating joined matrixtable with cohort info:")
    # annotate with cohorts
    mt_annotated = annotate_samples_with_cohort_info(mt, table_cohort)

    mt_annotated = mt_annotated.key_rows_by('locus').distinct_by_row().key_rows_by('locus', 'alleles')
    # mt_split = hl.split_multi_hts(
    #    mt, keep_star=False, left_aligned=False, permit_shuffle=True)
    # mt_split = mt_split.checkpoint(
    #   f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}-split-multi_cohorts.mt",  overwrite=True)
    print("Now writing joined matrixtable to disk:")
    mt_annotated.write(
        f"{lustre_dir}/MegaWESSanger_cohorts.mt", overwrite=True)
    print(f"Wrote matrixtable for whole genome.")
