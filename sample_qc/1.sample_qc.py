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
    CHROMOSOME = "WGS"
    mt = hl.read_matrix_table(
        f"{temp_dir}/ddd-elgh-ukbb/{CHROMOSOME}_annotated.mt")
    mt_split = hl.split_multi_hts(mt, keep_star=False, left_aligned=False)
    mt_split = mt_split.checkpoint(
        f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}-split-multi_checkpoint.mt",  overwrite=True)
    print("Finished splitting and writing mt. ")
    mt = mt_split.annotate_rows(
        Variant_Type=hl.cond((hl.is_snp(mt_split.alleles[0], mt_split.alleles[1])), "SNP",
                             hl.cond(
            hl.is_insertion(
                mt_split.alleles[0], mt_split.alleles[1]),
            "INDEL",
            hl.cond(hl.is_deletion(mt_split.alleles[0],
                                   mt_split.alleles[1]), "INDEL",
                    "Other"))))
    mt_sampleqc = hl.sample_qc(mt, name='sample_QC_Hail')
    panda_df_unfiltered_table = mt_sampleqc.cols().flatten()
    print("Sex imputation:")
    #mt2_sex = mt2.select_entries(GT=hl.unphased_diploid_gt_index_call(mt2.GT.n_alt_alleles()))
    imputed_sex = hl.impute_sex(mt_sampleqc.GT)

    # Annotate samples male or female:
    mt = mt_sampleqc.annotate_cols(sex=hl.cond(
        imputed_sex[mt_sampleqc.s].is_female, "female", "male"))

    print("Outputting table of sample qc")
    panda_df_unfiltered_table.export(
        f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}_sampleQC_unfiltered_sex_annotated.tsv.bgz", header=True)

   # mt2 = hl.variant_qc(mt_sampleqc, name='variant_QC_Hail')

    #print('Exporting variant qc pandas table to disk')
   # mt_rows = mt2.rows()
   # mt_rows.select(mt_rows.variant_QC_Hail).flatten().export(f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}_variantQC_unfiltered.tsv.bgz",
    #          header=True)

    # run sex determination script -not without chrX!

    mt = mt.checkpoint(
        f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}-sampleqc-unfiltered_sex_annotated.mt", overwrite=True)
    #
    # run sample_qc
    # plot various sample_qc per cohort -use intervalwgs threshold
    # plot sex determination per cohort - also undetermined -not without chrX
    # run variant_qc
