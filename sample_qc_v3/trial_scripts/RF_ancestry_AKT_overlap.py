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


def get_reference_genome(
    locus: Union[hl.expr.LocusExpression, hl.expr.IntervalExpression],
    add_sequence: bool = False,
) -> hl.ReferenceGenome:
    """
    Returns the reference genome associated with the input Locus expression
    :param locus: Input locus
    :param add_sequence: If set, the fasta sequence is added to the reference genome
    :return: Reference genome
    """
    if isinstance(locus, hl.expr.LocusExpression):
        ref = locus.dtype.reference_genome
    else:
        assert isinstance(locus, hl.expr.IntervalExpression)
        ref = locus.dtype.point_type.reference_genome
    if add_sequence:
        ref = add_reference_sequence(ref)
    return ref


def filter_to_autosomes(
    t: Union[hl.MatrixTable, hl.Table]
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filters the Table or MatrixTable to autosomes only.
    This assumes that the input contains a field named `locus` of type Locus
    :param t: Input MT/HT
    :return:  MT/HT autosomes
    """
    reference = get_reference_genome(t.locus)
    autosomes = hl.parse_locus_interval(
        f"{reference.contigs[0]}-{reference.contigs[21]}", reference_genome=reference
    )
    return hl.filter_intervals(t, [autosomes])


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
        "s3a://DDD-ELGH-UKBB-exomes/ancestry/sanger_cohort_known_populations_ukbb.tsv", delimiter="\t").key_by('s')

    # drop cohorts
    # annotate with cohorts and populations from s3 table.
    # save matrixtable
    mt = hl.read_matrix_table(
        f"{temp_dir}/ddd-elgh-ukbb/relatedness_ancestry/ddd-elgh-ukbb/chr1_chr20_XY_ldpruned.mt")
    mt = mt.annotate_cols(cohort=cohorts_pop[mt.s].cohort)
    mt = mt.annotate_cols(known_pop=cohorts_pop[mt.s].known_population)
    mt = mt.annotate_cols(gVCF=cohorts_pop[mt.s].gVCF_ID)
    # done the above on pca_RF jupyter notebook
    # mt = hl.read_matrix_table(
    #    f"{temp_dir}/ddd-elgh-ukbb/Sanger_cohorts_chr1-20-XY_new_cohorts.mt")
    #mt = hl.split_multi_hts(    mt, keep_star=False, left_aligned=False)
    mt.write(
        f"{tmp_dir}/ddd-elgh-ukbb/Sanger_chr1-20-XY_new_cohorts_split_multi_ld_pruned.mt", overwrite=True)
    # filter matrixtable
    logger.info("wrote mt ")
    # filter mt
    # mt = hl.read_matrix_table(
    #    f"{temp_dir}/ddd-elgh-ukbb/Sanger_cohorts_chr1-20-XY_new_cohorts_split_multi.mt")
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
    mt = mt.filter_rows(~ hl.is_mnp(mt.alleles[0], mt.alleles[1]))
    mt = mt.filter_rows(~ hl.is_indel(mt.alleles[0], mt.alleles[1]))
    mt = mt.filter_rows(~ hl.is_complex(mt.alleles[0], mt.alleles[1]))
    mt_vqc = hl.variant_qc(mt, name='variant_QC_Hail')
    # (mt_vqc.variant_QC_Hail.p_value_hwe >= 10 ** -6) & not to use this according to hcm.
    mt_vqc_filtered = mt_vqc.filter_rows(
        (mt_vqc.variant_QC_Hail.call_rate >= 0.99) &
        (mt_vqc.variant_QC_Hail.AF[1] >= 0.05) &
        (mt_vqc.variant_QC_Hail.AF[1] <= 0.95)
    )
    mt_vqc_filtered = mt_vqc_filtered.filter_rows(hl.is_defined(
        bed_to_exclude_pca[mt_vqc_filtered.locus]), keep=False)
    # overlap AKT dataset:

    mt_1kg_chr1_chr20 = hl.read_matrix_table(
        f"{temp_dir}/ddd-elgh-ukbb/relatedness_ancestry/ancestry_work/1000g_chr1_20_AKT_overlap.mt")
    #mt_vqc_filtered1 = mt_vqc_filtered.key_rows_by("locus")
    mt_1kg_chr1_chr20 = mt_1kg_chr1_chr20.key_rows_by("locus")
    mt_vqc_filtered = mt_vqc_filtered.filter_rows(
        hl.is_defined(mt_1kg_chr1_chr20.rows()[mt_vqc_filtered.locus]))
    logger.info("done filtering writing mt")
    # ld pruning
    pruned_ht = hl.ld_prune(mt_vqc_filtered.GT, r2=0.2, bp_window_size=500000)
    #pruned_ht = hl.ld_prune(mt.GT, r2=0.1)
    pruned_mt = mt_vqc_filtered.filter_rows(
        hl.is_defined(pruned_ht[mt_vqc_filtered.row_key]))
    # remove pruned areas that need to be removed

    pruned_mt.write(
        f"{tmp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_ldpruned_updated.mt", overwrite=True)
    # pruned_mt = hl.read_matrix_table(
    #    f"{temp_dir}/ddd-elgh-ukbb/relatedness_ancestry/ddd-elgh-ukbb/chr1_chr20_XY_ldpruned.mt")

    # related_samples_to_drop = hl.read_table(
    #    f"{temp_dir}/ddd-elgh-ukbb/relatedness_ancestry/ddd-elgh-ukbb/chr1_chr20_XY_related_samples_to_remove.ht")

    logger.info("run_pca_with_relateds")
    # pca_evals, pca_scores, pca_loadings = run_pca_with_relateds(
    #    pruned_mt, related_samples_to_drop, autosomes_only=True)
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(
        pruned_mt.GT, k=10, compute_loadings=True)
    pca_scores = pca_scores.annotate(
        known_pop=pruned_mt.cols()[pca_scores.s].known_pop)
    # mt = mt.annotate_cols(
    #    loadings=pca_loadings[mt_vqc_filtered.col_key].loadings)
    # mt = mt.annotate_cols(known_pop="unk")
    # pca_scores = pca_scores.annotate(known_pop="unk")
    pca_scores.write(
        f"{tmp_dir}/ddd-elgh-ukbb/pca_scores_after_pruning.ht", overwrite=True)
    pca_loadings.write(
        f"{tmp_dir}/ddd-elgh-ukbb/pca_loadings_after_pruning.ht", overwrite=True)
    with open(f"{temp_dir}/ddd-elgh-ukbb/pca_evals_after_pruning.txt", 'w') as f:
        for val in pca_evals:
            f.write(str(val))
    # pca_scores = hl.read_table(f"{temp_dir}/ddd-elgh-ukbb/pca_scores.ht")
    # pca_loadings = hl.read_table(f"{temp_dir}/ddd-elgh-ukbb/pca_loadings.ht")
    logger.info("assign population pcs")
   # population_assignment_table = assign_population_pcs(
    #    pca_scores, pca_loadings, known_col="known_pop")

    pop_ht, pop_clf = assign_population_pcs(
        pca_scores, pca_scores.scores, known_col="known_pop", n_estimators=100, prop_train=0.8)
    pop_ht.write(
        f"{tmp_dir}/ddd-elgh-ukbb/pop_assignments_test.ht", overwrite=True)
    pop_ht.export(f"{temp_dir}/ddd-elgh-ukbb/pop_assignments_test.tsv.gz")
    filename = f"{temp_dir}/ddd-elgh-ukbb/RF_model.pkl"
    pickle.dump(pop_clf, open(filename, 'wb'))
