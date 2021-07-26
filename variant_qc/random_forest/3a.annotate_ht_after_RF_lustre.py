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

project_root = Path(__file__).parent.parent.parent
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


tmp_dir = "hdfs://spark-master:9820/"
temp_dir = "file:///home/ubuntu/data/tmp"
plot_dir = "/home/ubuntu/data/tmp"
lustre_dir = "file:///lustre/scratch123/teams/hgi/mercury/megaWES-variantqc"


def count_trans_untransmitted_singletons(mt_filtered, ht):
    
    mt_trans = mt_filtered.filter_entries(mt_filtered.info.AC[0] == 2)
    mt_untrans = mt_filtered.filter_entries(mt_filtered.info.AC[0] == 1)
    
    mt_trans_count=mt_trans.group_cols_by(mt_trans.id).aggregate(transmitted_singletons_count=hl.agg.count_where(
                               # (mt_trans.info.AC[0] == 2) &
                                (mt_trans.proband_entry.GT.is_non_ref()) &
                                (
                                (mt_trans.father_entry.GT.is_non_ref()  )
                                 |
                                (mt_trans.mother_entry.GT.is_non_ref())
                                )
                                ))
    

    
    Total_transmitted_singletons=mt_trans_count.aggregate_entries(hl.agg.count_where(mt_trans_count.transmitted_singletons_count >0))
    print(Total_transmitted_singletons)
    mt_untrans_count = (mt_untrans.group_cols_by(mt_untrans.id).aggregate(
    untransmitted_singletons_count=hl.agg.count_where(
                   # (mt_untrans.info.AC[0] == 1) &
                    (mt_untrans.proband_entry.GT.is_hom_ref()) &
                    (
                    (mt_untrans.father_entry.GT.is_non_ref()) 
                    |
                    (mt_untrans.mother_entry.GT.is_non_ref())
                    )
                     )))
    Total_untransmitted_singletons=mt_untrans_count.aggregate_entries(hl.agg.count_where(mt_untrans_count.untransmitted_singletons_count >0))
    print(Total_untransmitted_singletons)
    Ratio_transmitted_untransmitted=Total_transmitted_singletons/Total_untransmitted_singletons
    print(Ratio_transmitted_untransmitted)

    mt2=mt_trans_count.annotate_rows(variant_transmitted_singletons=hl.agg.count_where(mt_trans_count.transmitted_singletons_count==1))
    mt2.variant_transmitted_singletons.summarize()

    mt3=mt_untrans_count.annotate_rows(variant_untransmitted_singletons=hl.agg.count_where(mt_untrans_count.untransmitted_singletons_count==1))
    mt3.variant_untransmitted_singletons.summarize()

    ht=ht.annotate(variant_transmitted_singletons=mt2.rows()[ht.key].variant_transmitted_singletons)
    ht=ht.annotate(variant_untransmitted_singletons=mt3.rows()[ht.key].variant_untransmitted_singletons)
    return(ht)
    


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

    #################################

   
    run_hash="91ba5f38"

    #ht = hl.read_table(
    #    f'{lustre_dir}/variant_qc/models/{run_hash}/rf_result_MegaWES_new.ht')
    ht = hl.read_table(
        f'{lustre_dir}/variant_qc/models/{run_hash}/rf_result_MegaWES_new.ht')
        
    #synonymous="s3a://DDD-ELGH-UKBB-exomes/vep/vep_consequences_worst_synonymous.tsv"
    #Look at notebooks/lustre_prepare_ht_after.RF.ipynb about how the following ht was created
    #Annotate with synonymous consequences
    ht_synonymous=hl.read_table( f'{lustre_dir}/grch38_synonymous_variants.ht')
    ht=ht.annotate(consequence=ht_synonymous[ht.key].consequence)
    print(ht.consequence.summarize())
     #annotate with de novo table
    ht_denovo=hl.read_table(
        f'{lustre_dir}/variant_qc/MegaWES_denovo_table.ht')
    ht=ht.annotate(de_novo_data=ht_denovo[ht.key].de_novo_data)

    #annotate with family stats
    ht_familystats=hl.read_table(f'{lustre_dir}/variant_qc/MegaWES_family_stats.ht')
    ht=ht.annotate(family_stats=ht_familystats[ht.key].family_stats)
    ht_stats=hl.read_table(f'{lustre_dir}/variant_qc/MegaWES_stats.ht')
    ht=ht.annotate(fam=ht_stats[ht.key])
    ht=ht.checkpoint(f'{lustre_dir}/variant_qc/models/{run_hash}_megaWES_RF_SYNONYMOUS_denovo_family_stats.ht', overwrite=True)
    
    
    #annotate with transmitted singletons counts
    mt_trios = hl.read_matrix_table(
        f'{lustre_dir}/variant_qc/MegaWES_trios_adj.mt')
    #mt_trios = hl.read_matrix_table(
    #    f'{lustre_dir}/MegaWES_trio_table.mt')
    mt_trios = mt_trios.annotate_rows(consequence=ht[mt_trios.row_key].consequence)



    mt_trios = mt_trios.checkpoint(
        f'{lustre_dir}/sanger_cohorts_trios_consequence.mt', overwrite=True)
    mt_filtered = mt_trios.filter_rows((mt_trios.info.AC[0] <= 2) & (
        mt_trios.consequence == "synonymous_variant"))
    #mt_filtered=hl.read_matrix_table(f'{lustre_dir}/variant_qc/MegaWESSanger_cohorts_AC_synonymous_filtered.mt')
    #mt_filtered=filter_to_autosomes(mt_filtered)
    mt_filtered=mt_filtered.checkpoint(f'{lustre_dir}/variant_qc/MegaWESSanger_cohorts_AC_synonymous_filtered_july_2021.mt',overwrite=True)
    
    ht=count_trans_untransmitted_singletons(mt_filtered, ht)
    
    
    #validated denovos
    
    ht_val_filtered=hl.read_table(f'{lustre_dir}/variant_qc/DDD_validated_denovo_b38_only_denovo_interitance.ht')
    ht=ht.annotate(validated_denovo_inheritance=ht_val_filtered[ht.key].inheritance)

    ht.write(f'{lustre_dir}/variant_qc/models/{run_hash}_rf_result_FINAL_for_RANKING.ht', overwrite=True)

    