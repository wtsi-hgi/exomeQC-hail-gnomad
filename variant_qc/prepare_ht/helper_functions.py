import hail as hl
from typing import Any, Counter, List, Optional, Tuple, Union

from gnomad.resources.grch38 import gnomad
from gnomad.utils.annotations import annotate_adj
from gnomad.utils.annotations import unphase_call_expr, add_variant_type
from gnomad.utils.annotations import annotate_adj, bi_allelic_expr, bi_allelic_site_inbreeding_expr
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.sample_qc.relatedness import (
    SIBLINGS,
    generate_sib_stats_expr,
    generate_trio_stats_expr,
)
from gnomad.utils.filtering import filter_to_adj

from hail import Table


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
    # mt = mt.filter_cols(mt.meta.high_quality)
    fam_ht = hl.import_fam(fam_file, delimiter="\t")
    mt = mt.annotate_cols(unrelated_sample=hl.is_missing(fam_ht[mt.s]))
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    mt = annotate_adj(mt)
    mt = mt.annotate_rows(
        ac_qc_samples_raw=hl.agg.sum(mt.GT.n_alt_alleles()),
        # ac_qc_samples_unrelated_raw=hl.agg.filter(~mt.meta.all_samples_related, hl.agg.sum(mt.GT.n_alt_alleles())),
        # ac_release_samples_raw=hl.agg.filter(mt.meta.release, hl.agg.sum(mt.GT.n_alt_alleles())),
        ac_qc_samples_adj=hl.agg.filter(
            mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())),
        # ac_qc_samples_unrelated_adj=hl.agg.filter(~mt.meta.all_samples_related & mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())),
        # ac_release_samples_adj=hl.agg.filter(mt.meta.release & mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())),
    )
    return mt.rows()


def get_truth_ht(omni, mills, thousand_genomes, hapmap) -> Table:
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
