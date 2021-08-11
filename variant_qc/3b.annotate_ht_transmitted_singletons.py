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


tmp_dir = "hdfs://spark-master:9820/"
temp_dir = "file:///home/ubuntu/data/tmp"
plot_dir = "/home/ubuntu/data/tmp"
lustre_dir = "file:///lustre/scratch123/teams/hgi/mercury/megaWES-variantqc"

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

if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory
  
    hl.init(sc=sc, tmp_dir=lustre_dir, local_tmpdir=lustre_dir, default_reference="GRCh38")

    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])

    ################################
    ### DO NOT RUN THIS, this is now part of 3a.annotate_ht_after_RF_lustre.py
    #################################

   
    run_hash = "ae281191"
    ht = hl.read_table(
        f'{lustre_dir}/variant_qc/models/{run_hash}_megaWES_RF_SYNONYMOUS_denovo_family_stats.ht')

    
    mt_filtered=hl.read_matrix_table(f'{lustre_dir}/variant_qc/MegaWESSanger_cohorts_AC_synonymous_filtered.mt')
    mt_filtered=filter_to_autosomes(mt_filtered)

    mt_trans = mt_filtered.filter_entries(mt_filtered.info.AC[0] == 2)
    mt_untrans = mt_filtered.filter_entries(mt_filtered.info.AC[0] == 1)
    
    mt_trans_count=mt_trans.group_cols_by(mt_trans.id).aggregate(transmitted_singletons_count=hl.agg.count_where((mt_trans.info.AC[0] == 2) &
                                (mt_trans.proband_entry.GT.is_het_ref()) &
                                (mt_trans.father_entry.GT.is_het_ref()) |
                                (mt_trans.mother_entry.GT.is_het_ref())))
    
    Total_transmitted_singletons=mt_trans_count.aggregate_entries(hl.agg.count_where(mt_trans_count.transmitted_singletons_count==1))
    print(Total_transmitted_singletons)
    mt_untrans_count = (mt_untrans.group_cols_by(mt_untrans.id).aggregate(
    untransmitted_singletons_count=hl.agg.count_where((mt_untrans.info.AC[0] == 1) &
                     (mt_untrans.proband_entry.GT.is_hom_ref()) &
                     (mt_untrans.father_entry.GT.is_het_ref()) |
                     (mt_untrans.mother_entry.GT.is_het_ref()))))
    Total_untransmitted_singletons=mt_untrans_count.aggregate_entries(hl.agg.count_where(mt_untrans_count.untransmitted_singletons_count==1))
    print(Total_untransmitted_singletons)
    Ratio_transmitted_untransmitted=Total_transmitted_singletons/Total_untransmitted_singletons
    print(Ratio_transmitted_untransmitted)

    mt2=mt_trans_count.annotate_rows(variant_transmitted_singletons=hl.agg.count_where(mt_trans_count.transmitted_singletons_count==1))
    mt2.variant_transmitted_singletons.summarize()

    mt3=mt_untrans_count.annotate_rows(variant_untransmitted_singletons=hl.agg.count_where(mt_untrans_count.untransmitted_singletons_count==1))
    mt3.variant_untransmitted_singletons.summarize()

    ht=ht.annotate(variant_transmitted_singletons=mt2.rows()[ht.key].variant_transmitted_singletons)
    ht=ht.annotate(variant_untransmitted_singletons=mt3.rows()[ht.key].variant_untransmitted_singletons)
    ht.write(f'{lustre_dir}/variant_qc/models/{run_hash}_rf_result_FINAL_for_RANKING.ht', overwrite=True)
    #ht.write(f'{lustre_dir}/variant_qc/models/{run_hash}_rf_result_transmitted_singletons_final.ht', overwrite=True)