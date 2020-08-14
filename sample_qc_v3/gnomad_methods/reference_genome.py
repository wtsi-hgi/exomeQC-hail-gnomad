import logging
from typing import List, Optional, Union

import hail as hl

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
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
