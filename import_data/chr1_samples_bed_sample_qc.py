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
intersection_bed = "s3a://DDD-ELGH-UKBB-exomes/agilentv5_UKB_intersection.bed"
union_bed = "s3a://DDD-ELGH-UKBB-exomes/agilentv5_UKB_union.bed"

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

    #####################################################################
    ###################### INPUT DATA  ##############################
    #####################################################################
    CHROMOSOME = "chr1"
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

    intersection_table = hl.import_bed(
        intersection_bed, reference_genome='GRCh38')
    union_table = hl.import_bed(union_bed, reference_genome='GRCh38')
    mt_intersection = mt.filter_rows(
        hl.is_defined(intersection_table[mt.locus]))
    mt_union = mt.filter_rows(hl.is_defined(union_table[mt.locus]))

    mt_intersection = hl.sample_qc(mt_intersection, name='sample_QC_Hail')
    pandadf1 = mt_intersection.cols().flatten()
    print("Outputting table of sample qc")
    pandadf1.export(
        f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}_intersection_BED_sampleQC.tsv.bgz", header=True)

    mt_intersection = mt_intersection.checkpoint(
        f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}-intersection_BED.mt", overwrite=True)

    mt = mt.checkpoint(
        f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}-sampleqc-unfiltered_sex_annotated.mt", overwrite=True)

    mt_union = hl.sample_qc(mt_union, name='sample_QC_Hail')
    pandadf2 = mt_union.cols().flatten()
    print("Outputting table of sample qc")
    pandadf2.export(
        f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}_union_BED_sampleQC.tsv.bgz", header=True)
    mt_union = mt_union.checkpoint(
        f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}-union_BED.mt", overwrite=True)
    #
    # run sample_qc
    # plot various sample_qc per cohort -use intervalwgs threshold
    # plot sex determination per cohort - also undetermined -not without chrX
    # run variant_qc
