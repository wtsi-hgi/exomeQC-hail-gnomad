# Pavlos Antoniou
# 16/09/2020
#  trio matrixtable creation from fam file
import os
import hail as hl
import pandas as pd
import pyspark
import json
import sys
import re
from pathlib import Path
import logging
from typing import Any, Counter, List, Optional, Tuple, Union, Dict
from bokeh.plotting import output_file, save, show
from gnomad.resources.grch38 import gnomad
from gnomad.utils.annotations import annotate_adj
from gnomad.utils.annotations import unphase_call_expr, add_variant_type
from gnomad.utils.annotations import annotate_adj, bi_allelic_expr, bi_allelic_site_inbreeding_expr
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.filtering import filter_to_adj
from gnomad.sample_qc.relatedness import (
    SIBLINGS,
    generate_sib_stats_expr,
    generate_trio_stats_expr,
)

from hail import Table

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

LABEL_COL = "rf_label"
TRAIN_COL = "rf_train"
PREDICTION_COL = "rf_prediction"
INFO_FEATURES = [
    "AS_QD",
    "AS_ReadPosRankSum",
    "AS_MQRankSum",
    "AS_SOR",
]  # Note: AS_SOR is currently in VQSR HT and named SOR in the VQSR split HT
FEATURES = [
    "InbreedingCoeff",
    "variant_type",
    "allele_type",
    "n_alt_alleles",
    "was_mixed",
    "has_star",
    "AS_QD",
    "AS_pab_max",
    "AS_MQRankSum",
    "AS_SOR",
    "AS_ReadPosRankSum",
]
TRUTH_DATA = ["hapmap", "omni", "mills", "kgp_phase1_hc"]
INBREEDING_COEFF_HARD_CUTOFF = -0.3


def annotate_freq(
    mt: hl.MatrixTable,
    sex_expr: Optional[hl.expr.StringExpression] = None,
    pop_expr: Optional[hl.expr.StringExpression] = None,
    subpop_expr: Optional[hl.expr.StringExpression] = None,
    additional_strata_expr: Optional[Dict[str, hl.expr.StringExpression]] = None,
    downsamplings: Optional[List[int]] = None,
) -> hl.MatrixTable:
    """
    Adds a row annotation `freq` to the input `mt` with stratified allele frequencies,
    and a global annotation `freq_meta` with metadata.

    .. note::

        Currently this only supports bi-allelic sites.
        The input `mt` needs to have the following entry fields:
        - GT: a CallExpression containing the genotype
        - adj: a BooleanExpression containing whether the genotype is of high quality or not.
        All expressions arguments need to be expression on the input `mt`.

    .. rubric:: `freq` row annotation

    The `freq` row annotation is an Array of Struct, with each Struct containing the following fields:

        - AC: int32
        - AF: float64
        - AN: int32
        - homozygote_count: int32

    Each element of the array corresponds to a stratification of the data, and the metadata about these annotations is
    stored in the globals.

    .. rubric:: Global `freq_meta` metadata annotation

    The global annotation `freq_meta` is added to the input `mt`. It is a list of dict.
    Each element of the list contains metadata on a frequency stratification and the index in the list corresponds
    to the index of that frequency stratification in the `freq` row annotation.

    .. rubric:: The `downsamplings` parameter

    If the `downsamplings` parameter is used, frequencies will be computed for all samples and by population
    (if `pop_expr` is specified) by downsampling the number of samples without replacement to each of the numbers specified in the
    `downsamplings` array, provided that there are enough samples in the dataset.
    In addition, if `pop_expr` is specified, a downsampling to each of the exact number of samples present in each population is added.
    Note that samples are randomly sampled only once, meaning that the lower downsamplings are subsets of the higher ones.

    :param mt: Input MatrixTable
    :param sex_expr: When specified, frequencies are stratified by sex. If `pop_expr` is also specified, then a pop/sex stratifiction is added.
    :param pop_expr: When specified, frequencies are stratified by population. If `sex_expr` is also specified, then a pop/sex stratifiction is added.
    :param subpop_expr: When specified, frequencies are stratified by sub-continental population. Note that `pop_expr` is required as well when using this option.
    :param additional_strata_expr: When specified, frequencies are stratified by the given additional strata found in the dict. This can e.g. be used to stratify by platform.
    :param downsamplings: When specified, frequencies are computed by downsampling the data to the number of samples given in the list. Note that if `pop_expr` is specified, downsamplings by population is also computed.
    :return: MatrixTable with `freq` annotation
    """

    if subpop_expr is not None and pop_expr is None:
        raise NotImplementedError(
            "annotate_freq requires pop_expr when using subpop_expr"
        )

    if additional_strata_expr is None:
        additional_strata_expr = {}

    _freq_meta_expr = hl.struct(**additional_strata_expr)
    if sex_expr is not None:
        _freq_meta_expr = _freq_meta_expr.annotate(sex=sex_expr)
    if pop_expr is not None:
        _freq_meta_expr = _freq_meta_expr.annotate(pop=pop_expr)
    if subpop_expr is not None:
        _freq_meta_expr = _freq_meta_expr.annotate(subpop=subpop_expr)

    # Annotate cols with provided cuts
    mt = mt.annotate_cols(_freq_meta=_freq_meta_expr)

    # Get counters for sex, pop and subpop if set
    cut_dict = {
        cut: hl.agg.filter(
            hl.is_defined(mt._freq_meta[cut]), hl.agg.counter(
                mt._freq_meta[cut])
        )
        for cut in mt._freq_meta
        if cut != "subpop"
    }
    if "subpop" in mt._freq_meta:
        cut_dict["subpop"] = hl.agg.filter(
            hl.is_defined(mt._freq_meta.pop) & hl.is_defined(
                mt._freq_meta.subpop),
            hl.agg.counter(
                hl.struct(subpop=mt._freq_meta.subpop, pop=mt._freq_meta.pop)
            ),
        )

    cut_data = mt.aggregate_cols(hl.struct(**cut_dict))
    sample_group_filters = []

    # Create downsamplings if needed
    if downsamplings is not None:
        # Add exact pop size downsampling if pops were provided
        if cut_data.get("pop"):
            downsamplings = list(
                set(downsamplings + list(cut_data.get("pop").values()))
            )  # Add the pops values if not in yet
            downsamplings = sorted(
                [x for x in downsamplings if x <= sum(
                    cut_data.get("pop").values())]
            )
        logger.info(
            f"Found {len(downsamplings)} downsamplings: {downsamplings}")

        # Shuffle the samples, then create a global index for downsampling
        # And a pop-index if pops were provided
        downsampling_ht = mt.cols()
        downsampling_ht = downsampling_ht.annotate(r=hl.rand_unif(0, 1))
        downsampling_ht = downsampling_ht.order_by(downsampling_ht.r)
        scan_expr = {"global_idx": hl.scan.count()}
        if cut_data.get("pop"):
            scan_expr["pop_idx"] = hl.scan.counter(downsampling_ht._freq_meta.pop).get(
                downsampling_ht._freq_meta.pop, 0
            )
        downsampling_ht = downsampling_ht.annotate(**scan_expr)
        downsampling_ht = downsampling_ht.key_by("s").select(*scan_expr)
        mt = mt.annotate_cols(downsampling=downsampling_ht[mt.s])
        mt = mt.annotate_globals(downsamplings=downsamplings)

        # Create downsampled sample groups
        sample_group_filters.extend(
            [
                (
                    {"downsampling": str(ds), "pop": "global"},
                    mt.downsampling.global_idx < ds,
                )
                for ds in downsamplings
            ]
        )
        if cut_data.get("pop"):
            sample_group_filters.extend(
                [
                    (
                        {"downsampling": str(ds), "pop": pop},
                        (mt.downsampling.pop_idx < ds) & (
                            mt._freq_meta.pop == pop),
                    )
                    for ds in downsamplings
                    for pop, pop_count in cut_data.get("pop", {}).items()
                    if ds <= pop_count
                ]
            )

    # Add all desired strata, starting with the full set and ending with downsamplings (if any)
    sample_group_filters = (
        [({}, True)]
        + [({"pop": pop}, mt._freq_meta.pop == pop)
            for pop in cut_data.get("pop", {})]
        + [({"sex": sex}, mt._freq_meta.sex == sex)
            for sex in cut_data.get("sex", {})]
        + [
            (
                {"pop": pop, "sex": sex},
                (mt._freq_meta.sex == sex) & (mt._freq_meta.pop == pop),
            )
            for sex in cut_data.get("sex", {})
            for pop in cut_data.get("pop", {})
        ]
        + [
            (
                {"subpop": subpop.subpop, "pop": subpop.pop},
                (mt._freq_meta.pop == subpop.pop)
                & (mt._freq_meta.subpop == subpop.subpop),
            )
            for subpop in cut_data.get("subpop", {})
        ]
        + [
            ({strata: str(s_value)}, mt._freq_meta[strata] == s_value)
            for strata in additional_strata_expr
            for s_value in cut_data.get(strata, {})
        ]
        + sample_group_filters
    )

    # Annotate columns with group_membership
    mt = mt.annotate_cols(group_membership=[x[1]
                                            for x in sample_group_filters])

    # Create and annotate global expression with meta information
    freq_meta_expr = [
        dict(**sample_group[0], group="adj") for sample_group in sample_group_filters
    ]
    freq_meta_expr.insert(1, {"group": "raw"})
    mt = mt.annotate_globals(freq_meta=freq_meta_expr)

    # Create frequency expression array from the sample groups
    freq_expr = hl.agg.array_agg(
        lambda i: hl.agg.filter(
            mt.group_membership[i] & mt.adj, hl.agg.call_stats(
                mt.GT, mt.alleles)
        ),
        hl.range(len(sample_group_filters)),
    )

    # Insert raw as the second element of the array
    freq_expr = (
        freq_expr[:1]
        .extend([hl.agg.call_stats(mt.GT, mt.alleles)])
        .extend(freq_expr[1:])
    )

    # Select non-ref allele (assumes bi-allelic)
    freq_expr = freq_expr.map(
        lambda cs: cs.annotate(
            AC=cs.AC[1],
            AF=cs.AF[
                1
            ],  # TODO This is NA in case AC and AN are 0 -- should we set it to 0?
            homozygote_count=cs.homozygote_count[1],
        )
    )

    # Return MT with freq row annotation
    return mt.annotate_rows(freq=freq_expr).drop("_freq_meta")


def get_truth_ht() -> Table:
    """
    Returns a table with the following annotations from the latest version of the corresponding truth data:
    - hapmap
    - kgp_omni (1000 Genomes intersection Onni 2.5M array)
    - kgp_phase_1_hc (high confidence sites in 1000 genonmes)
    - mills (Mills & Devine indels)
    :return: A table with the latest version of popular truth data annotations
    """
    omni_ht = hl.read_table(omni)
    mills_ht = hl.read_table(mills)
    thousand_genomes_ht = hl.read_table(thousand_genomes)
    hapmap_ht = hl.read_table(hapmap)
    return (
        hapmap_ht
        .select(hapmap=True)
        .join(omni_ht.select(omni=True), how="outer")
        .join(thousand_genomes_ht.select(kgp_phase1_hc=True), how="outer")
        .join(mills_ht.select(mills=True), how="outer")
        .repartition(200, shuffle=False)
        .persist()
    )


def generate_trio_stats(
    mt: hl.MatrixTable, autosomes_only: bool = True, bi_allelic_only: bool = True
) -> hl.Table:
    """
    Default function to run `generate_trio_stats_expr` to get trio stats stratified by raw and adj
    .. note::
        Expects that `mt` is it a trio matrix table that was annotated with adj and if dealing with
        a sparse MT `hl.experimental.densify` must be run first.
        By default this pipeline function will filter `mt` to only autosomes and bi-allelic sites.
    :param mt: A Trio Matrix Table returned from `hl.trio_matrix`. Must be dense
    :param autosomes_only: If set, only autosomal intervals are used.
    :param bi_allelic_only: If set, only bi-allelic sites are used for the computation
    :return: Table with trio stats
    """
    if autosomes_only:
        mt = filter_to_autosomes(mt)
    if bi_allelic_only:
        mt = mt.filter_rows(bi_allelic_expr(mt))

    logger.info(f"Generating trio stats using {mt.count_cols()} trios.")
    trio_adj = mt.proband_entry.adj & mt.father_entry.adj & mt.mother_entry.adj

    ht = mt.select_rows(
        **generate_trio_stats_expr(
            mt,
            transmitted_strata={"raw": True, "adj": trio_adj},
            de_novo_strata={"raw": True, "adj": trio_adj},
            ac_strata={"raw": True, "adj": trio_adj},
        )
    ).rows()

    return ht


def generate_allele_data(mt: hl.MatrixTable) -> hl.Table:
    """
    Writes bi-allelic sites MT with the following annotations:
     - allele_data (nonsplit_alleles, has_star, variant_type, and n_alt_alleles)

    :param MatrixTable mt: Full unsplit MT
    :return: Table with allele data annotations
    :rtype: Table
    """
    ht = mt.rows().select()
    allele_data = hl.struct(nonsplit_alleles=ht.alleles,
                            has_star=hl.any(lambda a: a == '*', ht.alleles))
    ht = ht.annotate(allele_data=allele_data.annotate(
        **add_variant_type(ht.alleles)))

    ht = hl.split_multi_hts(ht)
    allele_type = (hl.case()
                   .when(hl.is_snp(ht.alleles[0], ht.alleles[1]), 'snv')
                   .when(hl.is_insertion(ht.alleles[0], ht.alleles[1]), 'ins')
                   .when(hl.is_deletion(ht.alleles[0], ht.alleles[1]), 'del')
                   .default('complex')
                   )
    ht = ht.annotate(allele_data=ht.allele_data.annotate(allele_type=allele_type,
                                                         was_mixed=ht.allele_data.variant_type == 'mixed'))
    return ht


def generate_ac(mt: hl.MatrixTable, fam_file: str) -> hl.Table:
    """
    Creates Table with QC samples, QC samples removing children and release samples raw and adj ACs.
    """
    #mt = mt.filter_cols(mt.meta.high_quality)
    fam_ht = hl.import_fam(fam_file, delimiter="\t")
    mt = mt.annotate_cols(unrelated_sample=hl.is_missing(fam_ht[mt.s]))
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    mt = annotate_adj(mt)
    mt = mt.annotate_rows(
        ac_qc_samples_raw=hl.agg.sum(mt.GT.n_alt_alleles()),
        #ac_qc_samples_unrelated_raw=hl.agg.filter(~mt.meta.all_samples_related, hl.agg.sum(mt.GT.n_alt_alleles())),
        #ac_release_samples_raw=hl.agg.filter(mt.meta.release, hl.agg.sum(mt.GT.n_alt_alleles())),
        ac_qc_samples_adj=hl.agg.filter(
            mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())),
        #ac_qc_samples_unrelated_adj=hl.agg.filter(~mt.meta.all_samples_related & mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())),
        #ac_release_samples_adj=hl.agg.filter(mt.meta.release & mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())),
    )
    return mt.rows()


def generate_family_stats(mt: hl.MatrixTable, fam_file: str, calculate_adj: bool = False) -> Tuple[hl.Table, hl.Table]:
    """
    Writes bi-allelic sites MT with the following annotations:
     - family_stats (TDT, Mendel Errors, AC_unrelated_qc)
     - truth_data (presence in Omni, HapMap, 1KG high conf SNVs, Mills)

    :param MatrixTable mt: Full MT
    :param str fam_file: Fam pedigree file location
    :param bool calculate_adj: Whether to also calculate family metrics for adj genotypes
    :return: Table with qc annotations
    :rtype: Table
    """
    #mt = mt.select_cols(high_quality=mt.meta.high_quality)
    mt = mt.select_rows()
    mt = annotate_unrelated_sample(mt, fam_file)

    # Unphased for now, since mendel_errors does not support phased alleles
    mt = mt.annotate_entries(GT=unphase_call_expr(mt.GT))
    ped = hl.Pedigree.read(fam_file, delimiter='\\t')
    family_stats_struct, family_stats_sample_ht = family_stats(mt, ped, 'raw')
    mt = mt.annotate_rows(family_stats=[family_stats_struct])

    if calculate_adj:
        mt = filter_to_adj(mt)
        adj_family_stats_struct, adj_family_stats_sample_ht = family_stats(
            mt, ped, 'adj')

        family_stats_sample_ht = family_stats_sample_ht.annotate(
            adj=adj_family_stats_sample_ht[family_stats_sample_ht.s])

        mt = mt.annotate_rows(
            family_stats=mt.family_stats.append(adj_family_stats_struct))

    return mt.rows(), family_stats_sample_ht


def family_stats(mt: hl.MatrixTable, ped: hl.Pedigree, group_name: str) -> Tuple[hl.expr.StructExpression, hl.Table]:
    tdt_table = hl.transmission_disequilibrium_test(mt, ped)
    _, _, per_sample, per_variant = hl.mendel_errors(mt.GT, ped)
    family_stats_struct = hl.struct(mendel=per_variant[mt.row_key],
                                    tdt=tdt_table[mt.row_key],
                                    unrelated_qc_callstats=hl.agg.filter(mt.unrelated_sample,
                                                                         hl.agg.call_stats(mt.GT, mt.alleles)),
                                    meta={'group': group_name})
    return family_stats_struct, per_sample


def read_fam(fam_file: str) -> hl.Table:
    columns = ['fam_id', 's', 'pat_id', 'mat_id', 'is_female']
    return hl.import_table(fam_file, no_header=True).rename({f'f{i}': c for i, c in enumerate(columns)}).key_by('s')


def annotate_unrelated_sample(mt: hl.MatrixTable, fam_file: str) -> hl.MatrixTable:
    fam_ht = read_fam(fam_file)
    return mt.annotate_cols(unrelated_sample=hl.is_missing(fam_ht[mt.s]))


def generate_de_novos(mt: hl.MatrixTable, fam_file: str, freq_data: hl.Table) -> hl.Table:
    mt = mt.select_cols()
    fam_ht = read_fam(fam_file).key_by()
    fam_ht = fam_ht.select(
        s=[fam_ht.s, fam_ht.pat_id, fam_ht.mat_id]).explode('s').key_by('s')
    mt = mt.filter_cols(hl.is_defined(fam_ht[mt.s]))
    mt = mt.select_rows()
    mt = hl.split_multi_hts(mt)
    mt = mt.annotate_rows(family_stats=freq_data[mt.row_key].family_stats)
    ped = hl.Pedigree.read(fam_file, delimiter='\\t')

    de_novo_table = hl.de_novo(
        mt, ped, mt.family_stats[0].unrelated_qc_callstats.AF[1])
    de_novo_table = de_novo_table.key_by(
        'locus', 'alleles').collect_by_key('de_novo_data')

    return de_novo_table


if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory
    tmp_dir = "hdfs://spark-master:9820/"
    temp_dir = "file:///home/ubuntu/data/tmp"
    plot_dir = "/home/ubuntu/data/tmp"
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])

    ################################
    omni = f'{temp_dir}/ddd-elgh-ukbb/training_sets/1000G_omni2.5.hg38.ht'
    omni_ht = hl.read_table(omni)
    mills = f'{temp_dir}/ddd-elgh-ukbb/training_sets/Mills_and_1000G_gold_standard.indels.hg38.ht'
    mills_ht = hl.read_table(mills)
    thousand_genomes = f'{temp_dir}/ddd-elgh-ukbb/training_sets/1000G_phase1.snps.high_confidence.hg38.ht'
    thousand_genomes_ht = hl.read_table(thousand_genomes)
    hapmap = f'{temp_dir}/ddd-elgh-ukbb/training_sets/hapmap_3.3.hg38.ht'
    hapmap_ht = hl.read_table(hapmap)
    truthset_table = hl.read_table(
        f'{temp_dir}/ddd-elgh-ukbb/training_sets/truthset_table.ht')
    #################################

    # trio_stats_table = hl.read_table(
    #    f'{temp_dir}/ddd-elgh-ukbb/variant_qc/Sanger_cohorts_trios_stats.ht')
    group = "raw"

    mt = hl.read_matrix_table(
        f'{temp_dir}/ddd-elgh-ukbb/variant_qc/Sanger_cohorts_chr1-7and20_split.mt')
    fam = "s3a://DDD-ELGH-UKBB-exomes/trios/DDD_trios.fam"
    pedigree = hl.Pedigree.read(fam)
    (mt1, famstats_ht) = generate_family_stats(mt, fam)
    print("Writing mt and family stats_ht")
    mt1.write(f'{tmp_dir}/Sanger_cohorts_family_stats.mt', overwrite=True)
    famstats_ht.write(f'{tmp_dir}/Sanger_cohorts_family_stats.ht',overwrite=True)
    #mt = hl.variant_qc(mt)
    # mt=annotate_freq(mt)
    # trio annotation:
    #mt_adj = annotate_adj(mt)
    fam = "s3a://DDD-ELGH-UKBB-exomes/trios/DDD_trios.fam"
    pedigree = hl.Pedigree.read(fam)
    #trio_dataset = hl.trio_matrix(mt_adj, pedigree, complete_trios=True)
    trio_dataset = hl.read_matrix_table(
        f'{temp_dir}/ddd-elgh-ukbb/variant_qc/mt_trios_adj.mt')

    # trio_dataset=annotate_freq(trio_dataset)
    #trio_dataset = trio_dataset.annotate_rows(family_stats=freq_data[trio_dataset.row_key].family_stats)
    de_novo_table = hl.de_novo(
        mt, pedigree, trio_dataset.family_stats[0].unrelated_qc_callstats.AF[1])
    de_novo_table = de_novo_table.key_by(
        'locus', 'alleles').collect_by_key('de_novo_data')
    de_novo_table.write(
        f'{tmp_dir}/Sanger_cohorts_denovo_table.ht', overwrite=True)
