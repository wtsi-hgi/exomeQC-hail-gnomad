# Pavlos Antoniou
# 12/08/2020
#  QC for population assignment part 1
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
from typing import Any, Counter, List, Optional, Tuple, Union, Dict, Iterable
from sample_qc_v3.gnomad_methods.gnomad_filtering import compute_stratified_metrics_filter

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
tmp_dir = "hdfs://spark-master:9820/"
temp_dir = "file:///home/ubuntu/data/tmp"
plot_dir = "/home/ubuntu/data/tmp"
AKT_overlap = "s3a://DDD-ELGH-UKBB-exomes/ancestry/WES_AKT_1kg_intersection.vcf.mt"
cohorts_populations = "s3a://DDD-ELGH-UKBB-exomes/ancestry/sanger_cohort_known_populations_ukbb_elgh_labels_updated.tsv"
locations_exclude_from_pca = f"{temp_dir}/1000g/price_high_ld.bed.txt"

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


def main(args):

    # Read mt
    mt = hl.read_matrix_table(args.matrixtable)
    # pca_scores_pop
    pca_scores_pop = hl.read_table(args.pca_scores_population)

    # annotate mt with pop and superpop
    mt = mt.annotate_cols(assigned_pop=pca_scores_pop[mt.s].pop)

    # do sample_qc
    # calculate and annotate with metric heterozygosity
    mt_with_sampleqc = hl.sample_qc(mt, name='sample_qc')

    mt_with_sampleqc = mt_with_sampleqc.annotate_cols(sample_qc=mt_with_sampleqc.sample_qc.annotate(
        heterozygosity_rate=mt_with_sampleqc.sample_qc.n_het/mt_with_sampleqc.sample_qc.n_called))
    # save sample_qc and heterozygosity table as ht table
    mt_with_sampleqc.write(
        f"{args.output_dir}/ddd-elgh-ukbb/mt_pops_superpops_sampleqc.mt", overwrite=True)
    mt_with_sampleqc.cols().write(
        f"{args.output_dir}/ddd-elgh-ukbb/mt_pops_superpops_sampleqc.ht",  overwrite=True)
    pop_ht = hl.read_table(
        f"{args.output_dir}/ddd-elgh-ukbb/mt_pops_superpops_sampleqc.ht")
    # run function on metrics including heterozygosity first for pops:
    qc_metrics = ['heterozygosity_rate', 'n_snp', 'r_ti_tv',
                  'r_insertion_deletion', 'n_insertion', 'n_deletion', 'r_het_hom_var']
    pop_filter_ht = compute_stratified_metrics_filter(
        pop_ht, qc_metrics, ['assigned_pop'])
    pop_ht = pop_ht.annotate_globals(hl.eval(pop_filter_ht.globals))
    pop_ht = pop_ht.annotate(**pop_filter_ht[pop_ht.key]).persist()

    checkpoint = pop_ht.aggregate(hl.agg.count_where(
        hl.len(pop_ht.qc_metrics_filters) == 0))
    logger.info(f'{checkpoint} exome samples found passing pop filtering')
    pop_ht.write(f"{args.output_dir}/ddd-elgh-ukbb/mt_pops_QC_filters.ht")

    # run function on metrics including heterozygosity  for superpops:
    # pca_scores_superpop

    # pca_scores_superpop = hl.read_table(
    #    f"{temp_dir}/ddd-elgh-ukbb/new_labels/pop_assignments_updated_august2020_superpops.ht")
    #mt = mt.annotate_cols(assigned_superpop=pca_scores_superpop[mt.s].pop)
    # pop_ht_superpop = hl.read_table(
    #    f"{tmp_dir}/ddd-elgh-ukbb/mt_pops_superpops_sampleqc.ht")
    # pop_filter_ht = compute_stratified_metrics_filter(
    #    pop_ht_superpop, qc_metrics, ['assigned_superpop'])
    # pop_ht_superpop = pop_ht_superpop.annotate_globals(
    #    hl.eval(pop_filter_ht.globals))
    # pop_ht_superpop = pop_ht_superpop.annotate(
    #    **pop_filter_ht[pop_ht_superpop.key]).persist()

    # checkpoint = pop_ht_superpop.aggregate(hl.agg.count_where(
    #    hl.len(pop_ht_superpop.qc_metrics_filters) == 0))
    #logger.info(f'{checkpoint} exome samples found passing Superpop filtering')
    # pop_ht_superpop.write(
    #   f"{tmp_dir}/ddd-elgh-ukbb/mt_superpops_QC_filters.ht")


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

    #####################################################################
    ###################### INPUT DATA  ##############################
    #####################################################################
    parser = argparse.ArgumentParser()
    # Read the matrixtable, chrX and chrY should be included
    input_params = parser.add_argument_group("Input parameters")
    input_params.add_argument(
        "--matrixtable",
        help="Full path of input matrixtable after ld_pruning. Path format \"file:///home/ubuntu/data/tmp/path/to/.mt\"",
        default=f"{temp_dir}/ddd-elgh-ukbb/new_labels/chr1_chr20_ldpruned_updated.mt",
        type=str,
    )
    input_params.add_argument(
        "--pca_scores_population",
        help="Full path to hail table file .ht with PCA scores from previous PCA prediction script",
        default=f"{temp_dir}/ddd-elgh-ukbb/new_labels/pop_assignments_updated_august2020.ht",
        type=str
    )

    input_params.add_argument(
        "--output_dir",
        help="Full path of output folder to store results. Preferably hdfs or secure lustre",
        default=tmp_dir,
        type=str
    )

    args = parser.parse_args()
    main(args)
