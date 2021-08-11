
import logging
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import hail as hl

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def bi_allelic_site_inbreeding_expr(
    call: hl.expr.CallExpression,
) -> hl.expr.Float32Expression:
    """
    Return the site inbreeding coefficient as an expression to be computed on a MatrixTable.

    This is implemented based on the GATK InbreedingCoeff metric:
    https://software.broadinstitute.org/gatk/documentation/article.php?id=8032

    .. note::

        The computation is run based on the counts of alternate alleles and thus should only be run on bi-allelic sites.

    :param call: Expression giving the calls in the MT
    :return: Site inbreeding coefficient expression
    """

    def inbreeding_coeff(
        gt_counts: hl.expr.DictExpression,
    ) -> hl.expr.Float32Expression:
        n = gt_counts.get(0, 0) + gt_counts.get(1, 0) + gt_counts.get(2, 0)
        p = (2 * gt_counts.get(0, 0) + gt_counts.get(1, 0)) / (2 * n)
        q = (2 * gt_counts.get(2, 0) + gt_counts.get(1, 0)) / (2 * n)
        return 1 - (gt_counts.get(1, 0) / (2 * p * q * n))

    return hl.bind(inbreeding_coeff, hl.agg.counter(call.n_alt_alleles()))


def bi_allelic_expr(t: Union[hl.Table, hl.MatrixTable]) -> hl.expr.BooleanExpression:
    """
    Returns a boolean expression selecting bi-allelic sites only,
    accounting for whether the input MT/HT was split.

    :param t: Input HT/MT
    :return: Boolean expression selecting only bi-allelic sites
    """
    return ~t.was_split if "was_split" in t.row else (hl.len(t.alleles) == 2)
