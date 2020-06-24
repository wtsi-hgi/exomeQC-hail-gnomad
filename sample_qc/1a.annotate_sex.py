import os
import hail as hl
import pyspark
import json
import sys
import re
from pathlib import Path
import logging
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("unified_sample_qc_a")
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


def annotate_sex(mt: hl.MatrixTable, out_internal_mt_prefix: str,
                 male_threshold: float = 0.8, female_threshold: float = 0.5) -> hl.MatrixTable:
    """
    Imputes sex, exports data, and annotates mt with this data
    NOTE: Evaluated in R (plots) and decided on cutoff of F<0.5 for females and F>0.8 for males (default) for genomes

    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param str out_internal_mt_prefix: file path prefix for tsv containing samples and sex imputation annotations
    :return: MatrixTable with imputed sex annotations stashed in column annotation 'sex_check'
    :rtype: MatrixTable
    """
    mt1 = hl.filter_intervals(mt, [hl.parse_locus_interval('chrX')])
    #mt = mt.filter_rows(mt.locus.in_x_nonpar())
    mtx_unphased = mt1.select_entries(
        GT=hl.unphased_diploid_gt_index_call(mt1.GT.n_alt_alleles()))
    #imputed_sex = hl.impute_sex(mtx_unphased.GT)
    sex_ht = hl.impute_sex(mtx_unphased.GT, aaf_threshold=0.05,
                           female_threshold=female_threshold, male_threshold=male_threshold)
    sex_ht.export(out_internal_mt_prefix + '.sex_check.txt.bgz')
    sex_colnames = ['f_stat', 'is_female']
    sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    return mt


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
        f"{temp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_sex.mt")

   # From gnomad pply hard filters:


#    qc_ht = qc_ht.annotate(ambiguous_sex=((qc_ht.f_stat >= 0.5) & (hl.is_defined(qc_ht.normalized_y_coverage) & (qc_ht.normalized_y_coverage <= 0.1))) |
 #                          (hl.is_missing(qc_ht.f_stat)) |
  #                         ((qc_ht.f_stat >= 0.4) & (qc_ht.f_stat <= 0.6) & (hl.is_defined(
   #                            qc_ht.normalized_y_coverage) & (qc_ht.normalized_y_coverage > 0.1))),
    #                       sex_aneuploidy=(qc_ht.f_stat < 0.4) & hl.is_defined(qc_ht.normalized_y_coverage) & (qc_ht.normalized_y_coverage > 0.1))

  #  print("Annotating samples failing hard filters:")
   # logger.info("Annotating samples failing hard filters...")

    sex_expr = (hl.case()
                .when(mt.is_female, "female")
                .default("male"))

    mt = mt.annotate_cols(
        sex=sex_expr, data_type='exomes')
    mt.write(f"{tmp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_sex_annotations.mt",
             overwrite=True)

    # mt2_sex = mt2.select_entries(GT=hl.unphased_diploid_gt_index_call(mt2.GT.n_alt_alleles()))
    # imputed_sex = hl.impute_sex(mt_sampleqc.GT)

    # Annotate samples male or female:
    # mt = mt_sampleqc.annotate_cols(sex=hl.cond(
    #    imputed_sex[mt_sampleqc.s].is_female, "female", "male"))

    # print("Outputting table of sample qc")
    # panda_df_unfiltered_table.export(
    #    f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}_sampleQC_unfiltered_sex_annotated.tsv.bgz", header=True)

   # mt2 = hl.variant_qc(mt_sampleqc, name='variant_QC_Hail')

    # print('Exporting variant qc pandas table to disk')
   # mt_rows = mt2.rows()
   # mt_rows.select(mt_rows.variant_QC_Hail).flatten().export(f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}_variantQC_unfiltered.tsv.bgz",
    #          header=True)

    # run sex determination script -not without chrX!

    # mt = mt.checkpoint(
    #    f"{tmp_dir}/ddd-elgh-ukbb/{CHROMOSOME}-sampleqc-unfiltered_sex_annotated.mt", overwrite=True)
    #
    # run sample_qc
    # plot various sample_qc per cohort -use intervalwgs threshold
    # plot sex determination per cohort - also undetermined -not without chrX
    # run variant_qc
