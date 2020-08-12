import pandas as pd
import hail as hl
import logging
import random
from typing import Any, Counter, List, Optional, Tuple, Union, Dict, Iterable


def compute_stratified_metrics_filter(
    ht: hl.Table,
    qc_metrics: Dict[str, hl.expr.NumericExpression],
    strata: Optional[Dict[str, hl.expr.Expression]] = None,
    lower_threshold: float = 4.0,
    upper_threshold: float = 4.0,
    metric_threshold: Optional[Dict[str, Tuple[float, float]]] = None,
    filter_name: str = "qc_metrics_filters",
) -> hl.Table:
    """
    Compute median, MAD, and upper and lower thresholds for each metric used in outlier filtering

    :param ht: HT containing relevant sample QC metric annotations
    :param qc_metrics: list of metrics (name and expr) for which to compute the critical values for filtering outliers
    :param strata: List of annotations used for stratification. These metrics should be discrete types!
    :param lower_threshold: Lower MAD threshold
    :param upper_threshold: Upper MAD threshold
    :param metric_threshold: Can be used to specify different (lower, upper) thresholds for one or more metrics
    :param filter_name: Name of resulting filters annotation
    :return: Table grouped by strata, with upper and lower threshold values computed for each sample QC metric
    """

    _metric_threshold = {
        metric: (lower_threshold, upper_threshold) for metric in qc_metrics
    }
    if metric_threshold is not None:
        _metric_threshold.update(metric_threshold)

    def make_filters_expr(
        ht: hl.Table, qc_metrics: Iterable[str]
    ) -> hl.expr.SetExpression:
        return hl.set(
            hl.filter(
                lambda x: hl.is_defined(x),
                [hl.or_missing(ht[f"fail_{metric}"], metric)
                 for metric in qc_metrics],
            )
        )

    if strata is None:
        strata = {}

    ht = ht.select(**qc_metrics, **strata).key_by("s").persist()

    agg_expr = hl.struct(
        **{
            metric: hl.bind(
                lambda x: x.annotate(
                    lower=x.median - _metric_threshold[metric][0] * x.mad,
                    upper=x.median + _metric_threshold[metric][1] * x.mad,
                ),
                get_median_and_mad_expr(ht[metric]),
            )
            for metric in qc_metrics
        }
    )

    if strata:
        ht = ht.annotate_globals(
            qc_metrics_stats=ht.aggregate(
                hl.agg.group_by(hl.tuple([ht[x] for x in strata]), agg_expr),
                _localize=False,
            )
        )
        metrics_stats_expr = ht.qc_metrics_stats[hl.tuple(
            [ht[x] for x in strata])]
    else:
        ht = ht.annotate_globals(
            qc_metrics_stats=ht.aggregate(agg_expr, _localize=False)
        )
        metrics_stats_expr = ht.qc_metrics_stats

    fail_exprs = {
        f"fail_{metric}": (ht[metric] <= metrics_stats_expr[metric].lower)
        | (ht[metric] >= metrics_stats_expr[metric].upper)
        for metric in qc_metrics
    }
    ht = ht.transmute(**fail_exprs)
    stratified_filters = make_filters_expr(ht, qc_metrics)
    return ht.annotate(**{filter_name: stratified_filters})
