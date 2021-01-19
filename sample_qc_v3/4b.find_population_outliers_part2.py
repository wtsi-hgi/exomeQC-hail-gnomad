# Pavlos Antoniou
# 12/08/2020
#  QC for population assignment
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
from gnomad_filtering import compute_stratified_metrics_filter

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
        "s3a://DDD-ELGH-UKBB-exomes/ancestry/sanger_cohort_known_populations_ukbb_elgh_labels_updated.tsv", delimiter="\t").key_by('s')
    # Read mt
    mt = hl.read_matrix_table(
        f"{temp_dir}/ddd-elgh-ukbb/new_labels/chr1_chr20_ldpruned_updated.mt")
    # pca_scores_pop
    pca_scores_pop = hl.read_table(
        f"{temp_dir}/ddd-elgh-ukbb/new_labels/pop_assignments_updated_august2020.ht")

    strata = {}
    mt_with_sampleqc = hl.read_matrix_table(
        f"{tmp_dir}/ddd-elgh-ukbb/mt_pops_superpops_sampleqc.mt")
    pop_ht = hl.read_table(
        f"{tmp_dir}/ddd-elgh-ukbb/mt_pops_superpops_sampleqc.ht")
    # run function on metrics including heterozygosity first for pops:
    strata['assigned_pop'] = pop_ht.assigned_pop

    qc_metrics = {
        'sample_qc.heterozygosity_rate': pop_ht.sample_qc.heterozygosity_rate,
        'sample_qc.n_snp': pop_ht.sample_qc.n_snp,
        'sample_qc.r_ti_tv': pop_ht.sample_qc.r_ti_tv,
        'sample_qc.r_insertion_deletion': pop_ht.sample_qc.r_insertion_deletion,
        'sampleqc.n_insertion': pop_ht.sample_qc.n_insertion,
        'sampleqc.n_deletion': pop_ht.sample_qc.n_deletion,
        'sample_qc.r_het_hom_var': pop_ht.sample_qc.r_het_hom_var
    }
    pop_filter_ht = compute_stratified_metrics_filter(
        pop_ht, qc_metrics, strata)

    #pop_ht = pop_ht.annotate_globals(hl.eval(pop_filter_ht.globals))
    #pop_ht = pop_ht.annotate(**pop_filter_ht[pop_ht.key]).persist()

    # checkpoint = pop_ht.aggregate(hl.agg.count_where(
    #   hl.len(pop_ht.qc_metrics_filters) == 0))
    #logger.info(f'{checkpoint} exome samples found passing pop filtering')
    pop_filter_ht.write(
        f"{tmp_dir}/ddd-elgh-ukbb/mt_pops_QC_filters.ht", overwrite=True)

    # run function on metrics including heterozygosity  for superpops:
    strata = {}
    strata['assigned_superpop'] = pop_ht.assigned_superpop
    pop_filter_ht_superpop = compute_stratified_metrics_filter(
        pop_ht, qc_metrics, strata)

    pop_filter_ht_superpop.write(
        f"{tmp_dir}/ddd-elgh-ukbb/mt_superpops_QC_filters.ht", overwrite=True)
