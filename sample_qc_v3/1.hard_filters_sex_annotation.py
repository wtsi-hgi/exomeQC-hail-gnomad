# Pavlos Antoniou
# Sample QC pipeline first step
# Apply gnomad hard filters and calculate sex for each sample using hail's build in functions
# 19/01/2021
#  Required: a list of known population label for each sample

import os
import hail as hl
import pandas as pd
import pyspark
import argparse
import json
import sys
import re
from pathlib import Path
import logging
from typing import Any, Counter, List, Optional, Tuple, Union
from bokeh.plotting import output_file, save, show

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
tmp_dir = "hdfs://spark-master:9820/"
temp_dir = "file:///home/ubuntu/data/tmp"
plot_dir = "/home/ubuntu/data/tmp"


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


def main(args):
    mt = hl.read_matrix_table(args.matrixtable)

    # From gnomad apply hard filters:
    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                        (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > 0.001) &
                        (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99))
    mt.annotate_cols(callrate=hl.agg.fraction(hl.is_defined(mt.GT))).write(
        f"{args.output-dir}/mt_hard_filters_annotated.mt", overwrite=True)

    print("Sex imputation:")

    mt_sex = annotate_sex(
        mt, f"{args.output-dir}/sex_annotated", male_threshold=0.6)
    mt_sex.write(f"{args.output-dir}/mt_sex_annotated.mt", overwrite=True)

    qc_ht = mt_sex.cols()

    qc_ht = qc_ht.annotate(ambiguous_sex=((qc_ht.f_stat >= 0.5) & (hl.is_defined(qc_ht.normalized_y_coverage) & (qc_ht.normalized_y_coverage <= 0.1))) |
                           (hl.is_missing(qc_ht.f_stat)) |
                           ((qc_ht.f_stat >= 0.4) & (qc_ht.f_stat <= 0.6) & (hl.is_defined(
                               qc_ht.normalized_y_coverage) & (qc_ht.normalized_y_coverage > 0.1))),
                           sex_aneuploidy=(qc_ht.f_stat < 0.4) & hl.is_defined(qc_ht.normalized_y_coverage) & (qc_ht.normalized_y_coverage > 0.1))

    print("Annotating samples failing hard filters:")
    logger.info("Annotating samples failing hard filters...")

    sex_expr = (hl.case()
                .when(qc_ht.ambiguous_sex, "ambiguous_sex")
                .when(qc_ht.sex_aneuploidy, "sex_aneuploidy")
                .when(qc_ht.is_female, "female")
                .default("male"))

    qc_ht = qc_ht.annotate(
        sex=sex_expr, data_type='exomes').key_by('data_type', 's')
    qc_ht.write(f"{args.output-dir}/mt_ambiguous_sex_samples.ht",
                overwrite=True)


if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory

    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])
    parser = argparse.ArgumentParser()
    # Read the matrixtable, chrX and chrY should be included
    input_params = parser.add_argument_group("Input parameters")
    input_params.add_argument(
        "--matrixtable",
        help="Full path of input matrixtable. chrX and chrY variation should be included",
        default=f"{temp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_cohorts_split.mt",
        type=str,
    )
    input_params.add_argument(
        "--output-dir",
        help="Full path of output folder to store results. Preferably hdfs or secure lustre",
        default=tmp_dir,
        type=str
    )

    args = parser.parse_args()
    main(args)
