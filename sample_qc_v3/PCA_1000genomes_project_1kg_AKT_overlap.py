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


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


POP_NAMES = {
    "afr": "African/African-American",
    "ami": "Amish",
    "amr": "Latino",
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "eur": "European",
    "fin": "Finnish",
    "mde": "Middle Eastern",
    "nfe": "Non-Finnish European",
    "oth": "Other",
    "sas": "South Asian",
    "uniform": "Uniform",
    "sas_non_consang": "South Asian (F < 0.05)",
    "consanguineous": "South Asian (F > 0.05)",
    "exac": "ExAC",
    "bgr": "Bulgarian (Eastern European)",
    "est": "Estonian",
    "gbr": "British",
    "nwe": "North-Western European",
    "seu": "Southern European",
    "swe": "Swedish",
    "kor": "Korean",
    "sgp": "Singaporean",
    "jpn": "Japanese",
    "oea": "Other East Asian",
    "oeu": "Other European",
    "onf": "Other Non-Finnish European",
    "unk": "Unknown",
}

POP_COLORS = {
    "afr": "#941494",
    "ami": "#FFC0CB",
    "amr": "#ED1E24",
    "asj": "#FF7F50",
    "eas": "#108C44",
    "eur": "#6AA5CD",
    "fin": "#002F6C",
    "mde": "#33CC33",
    "nfe": "#6AA5CD",
    "oth": "#ABB9B9",
    "sas": "#FF9912",
    "uniform": "pink",
    "consanguineous": "pink",
    "sas_non_consang": "orange",
    "exac": "gray",
    "bgr": "#66C2A5",
    "est": "black",
    "gbr": "#C60C30",
    "nwe": "#C60C30",
    "seu": "#3CA021",
    "swe": "purple",
    "kor": "#4891D9",
    "sgp": "darkred",
    "jpn": "#BC002D",
    "oea": "#108C44",
    "oeu": "#6AA5CD",
    "onf": "#6AA5CD",
    "unk": "#ABB9B9",
    "": "#ABB9B9",
}


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


def pc_project(
    mt: hl.MatrixTable,
    loadings_ht: hl.Table,
    loading_location: str = "loadings",
    af_location: str = "pca_af",
) -> hl.Table:
    """
    Projects samples in `mt` on pre-computed PCs.

    :param mt: MT containing the samples to project
    :param loadings_ht: HT containing the PCA loadings and allele frequencies used for the PCA
    :param loading_location: Location of expression for loadings in `loadings_ht`
    :param af_location: Location of expression for allele frequency in `loadings_ht`
    :return: Table with scores calculated from loadings in column `scores`
    """
    n_variants = loadings_ht.count()

    mt = mt.annotate_rows(
        pca_loadings=loadings_ht[mt.row_key][loading_location],
        pca_af=loadings_ht[mt.row_key][af_location],
    )

    mt = mt.filter_rows(
        hl.is_defined(mt.pca_loadings)
        & hl.is_defined(mt.pca_af)
        & (mt.pca_af > 0)
        & (mt.pca_af < 1)
    )

    gt_norm = (mt.GT.n_alt_alleles() - 2 * mt.pca_af) / hl.sqrt(
        n_variants * 2 * mt.pca_af * (1 - mt.pca_af)
    )

    mt = mt.annotate_cols(scores=hl.agg.array_sum(mt.pca_loadings * gt_norm))

    return mt.cols().select("scores")


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
    project_mt = hl.read_matrix_table(
        f"{temp_dir}/ddd-elgh-ukbb/relatedness_ancestry/ddd-elgh-ukbb/chr1_chr20_XY_ldpruned.mt")
    vcf_1kg = "s3a://DDD-ELGH-UKBB-exomes/chr1_chr20_1kg_akt.vcf.gz"
    mt_1kg_chr1_chr20 = hl.import_vcf(vcf_1kg, array_elements_required=False,
                                      force_bgz=True)
    mt_1kg_chr1_chr20.write(
        f"{tmp_dir}/ddd-elgh-ukbb/1000g_chr1_20_AKT_overlap.mt", overwrite=True)

    # intersect this with the data we have
    logger.info("intersecting with project mt")
    project_mt = project_mt.key_rows_by("locus")
    mt_1kg = mt_1kg_chr1_chr20.filter_rows(
        hl.is_defined(project_mt.rows()[mt_1kg_chr1_chr20.locus]))
    mt_1kg = mt_1kg.checkpoint(
        f"{temp_dir}/ddd-elgh-ukbb/1000g_chr1_20_AKT_projrctdata_overlap.mt")
    #pruned_ht = hl.ld_prune(mt_vqc_filtered.GT, r2=0.2, bp_window_size=500000)

    # pruned_mt = hl.read_matrix_table(
    #    f"{temp_dir}/ddd-elgh-ukbb/1000g_chr1_20_snps_filtered_ldpruned.mt")
    # run pca
    logger.info("run pca")
    pca_evals, pca_scores, loadings_ht = hl.hwe_normalized_pca(
        mt_1kg.GT, k=10, compute_loadings=True)
    mt_1kg = mt_1kg.annotate_rows(
        af=hl.agg.mean(mt_1kg.GT.n_alt_alleles()) / 2)
    loadings_ht = loadings_ht.annotate(af=mt_1kg.rows()[loadings_ht.key].af)
    pca_scores.write(
        f"{tmp_dir}/ddd-elgh-ukbb/1K_AKT_data_pca_scores.ht", overwrite=True)
    loadings_ht.write(
        f"{tmp_dir}/ddd-elgh-ukbb/1K_AKT_data_pca_loadings.ht", overwrite=True)
    with open(f"{temp_dir}/ddd-elgh-ukbb/1K_AKT_data_pca_evals.txt", 'w') as f:
        for val in pca_evals:
            f.write(str(val))
    # pca_scores.export(f"{temp_dir}/ddd-elgh-ukbb/pca_scores.tsv.bgz")
    # loadings_ht.export(f"{temp_dir}/ddd-elgh-ukbb/pca_loadings.tsv.bgz")
       # ht = pc_project(
    #    project_mt, loadings_ht, "loadings", "af")

    ht = pc_project(
        project_mt, loadings_ht, "loadings", "af")

    ht.write(
        f"{tmp_dir}/ddd-elgh-ukbb/ourproject_1kg_AKT_scores.ht", overwrite=True)
    mt = project_mt.annotate_cols(scores=ht[project_mt.s].scores)
    mt.write(
        f"{tmp_dir}/ddd-elgh-ukbb/ourproject_1kg_AKT_projection.mt", overwrite=True)
    ht.export(f"{temp_dir}/ddd-elgh-ukbb/ourproject_1kg_AKT.txt.bgz")
