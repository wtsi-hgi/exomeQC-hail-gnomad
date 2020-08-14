import functools
import logging
import operator
from typing import Dict, List, Optional, Union


def filter_to_adj(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter genotypes to adj criteria
    """
    if "adj" not in list(mt.entry):
        mt = annotate_adj(mt)
    mt = mt.filter_entries(mt.adj)
    return mt.drop(mt.adj)


def filter_low_conf_regions(
    mt: Union[hl.MatrixTable, hl.Table],
    filter_lcr: bool = True,
    filter_decoy: bool = True,
    filter_segdup: bool = True,
    filter_exome_low_coverage_regions: bool = False,
    high_conf_regions: Optional[List[str]] = None,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filters low-confidence regions

    :param mt: MatrixTable or Table to filter
    :param filter_lcr: Whether to filter LCR regions
    :param filter_decoy: Whether to filter decoy regions
    :param filter_segdup: Whether to filter Segdup regions
    :param filter_exome_low_coverage_regions: Whether to filter exome low confidence regions
    :param high_conf_regions: Paths to set of high confidence regions to restrict to (union of regions)
    :return: MatrixTable or Table with low confidence regions removed
    """
    build = get_reference_genome(mt.locus).name
    if build == "GRCh37":
        import gnomad.resources.grch37.reference_data as resources
    elif build == "GRCh38":
        import gnomad.resources.grch38.reference_data as resources

    criteria = []
    if filter_lcr:
        lcr = resources.lcr_intervals.ht()
        criteria.append(hl.is_missing(lcr[mt.locus]))

    if filter_decoy:
        decoy = resources.decoy_intervals.ht()
        criteria.append(hl.is_missing(decoy[mt.locus]))

    if filter_segdup:
        segdup = resources.seg_dup_intervals.ht()
        criteria.append(hl.is_missing(segdup[mt.locus]))

    if filter_exome_low_coverage_regions:
        high_cov = resources.high_coverage_intervals.ht()
        criteria.append(hl.is_missing(high_cov[mt.locus]))

    if high_conf_regions is not None:
        for region in high_conf_regions:
            region = hl.import_locus_intervals(region)
            criteria.append(hl.is_defined(region[mt.locus]))

    if criteria:
        filter_criteria = functools.reduce(operator.iand, criteria)
        if isinstance(mt, hl.MatrixTable):
            mt = mt.filter_rows(filter_criteria)
        else:
            mt = mt.filter(filter_criteria)

    return mt


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
