# Pavlos Antoniou
# create Random Forest hail table
#  based on gnomad code
# 21/01/2021

import os
import hail as hl
import pandas as pd
import pyspark
import json
import sys
import re

from pathlib import Path
import logging
import argparse
from typing import Any, Counter, List, Optional, Tuple, Union
from bokeh.plotting import output_file, save, show
from gnomad.resources.grch38 import gnomad
from gnomad.utils.annotations import unphase_call_expr, add_variant_type
from gnomad.utils.reference_genome import add_reference_sequence
from gnomad.variant_qc.random_forest import (
    apply_rf_model,
    load_model,
    median_impute_features,
    pretty_print_runs,
    save_model,
)

from hail import Table

################################
# Define the hail persistent storage directory
tmp_dir = "hdfs://spark-master:9820/"
temp_dir = "file:///home/ubuntu/data/tmp"
plot_dir = "/home/ubuntu/data/tmp"
lustre_dir = "file:///lustre/scratch123/teams/hgi/mercury/megaWES-variantqc"
######################################
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
    "QD",
    "MQRankSum",
    "SOR",
    "ReadPosRankSum",
    "FS",
    "DP"
]
TRUTH_DATA = ["hapmap", "omni", "mills", "kgp_phase1_hc"]
INBREEDING_COEFF_HARD_CUTOFF = -0.3

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


def main(args):
    n_partitions = 500

    

    mt = hl.read_matrix_table(
        args.matrixtable)
    mt = mt.key_rows_by('locus').distinct_by_row(
    ).key_rows_by('locus', 'alleles')
   
    
    mt=filter_to_autosomes(mt)

    mt.write(f'{lustre_dir}/variant_qc/MegaWESSanger_cohorts_sampleQC_filtered_autosomes.mt', overwrite=True)


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

    #####################################################################
    ###################### INPUT DATA  ##############################
    #####################################################################
    parser = argparse.ArgumentParser()
    # Read the matrixtable
    input_params = parser.add_argument_group("Input parameters")
    input_params.add_argument(
        "--matrixtable",
        help="Full path of input matrixtable. Path format \"file:///home/ubuntu/data/tmp/path/to/.mt\"",
        default=f'{lustre_dir}/MegaWESSanger_cohorts_sampleQC_filtered.mt',
        type=str,
    )
    
    input_params.add_argument(
        "--output_dir",
        help="Full path of output folder to store results. Preferably hdfs or secure lustre",
        default=lustre_dir,
        type=str
    )
    args = parser.parse_args()
    main(args)
