
import os
import glob
from datetime import datetime
import hail as hl
import pyspark
import re
import pandas as pd
from pprint import pprint
from bokeh.plotting import output_file, save

from hail_methods.cutoffs import *
from hail_methods.hail_constants import *
from hail_methods.vep_hail_functions import *

class Timer(object):
    def __init__(self, message, now = datetime.utcnow):
        self.message = message
        self.now = now

    def __enter__(self):
        self.start = self.now()

    def __exit__(self, *exception):
        duration = self.now() - self.start

        status = "failed"
        if exception == (None, None, None):
            status = "completed"

        print(f"** {self.message} {status} after {duration}")


def split_mt_to_snps(mt:hl.MatrixTable) -> hl.MatrixTable:
    '''

    :param mt: hail matrixtable of all samples with both indels and SNVs
    :return: matrixtable with only the SNPs from all the samples
    '''

    mt_snps = hl.filter_alleles(mt, lambda allele, _: hl.is_snp(mt.alleles[0], allele))
    return mt_snps

def split_mt_to_indels(mt:hl.MatrixTable) -> hl.MatrixTable:
    '''

    :param mt: hail matrixtable of all samples with both indels and SNVs
    :return: hail matrixtable with only the indels
    '''
    mt_indels = hl.filter_alleles(mt,lambda allele, _: hl.is_indel(mt.alleles[0], allele))
    return mt_indels


def filter_snps_cutoffs(mt:hl.MatrixTable) -> hl.MatrixTable:
    '''
    Filter snps based on Hilary's cutoffs
    :param mt: snps only hail matrixtable
    :return: filtered matrixtable based on the  cutoffs listed in cutoffs.py file

    '''
    mt_snps_filtered = mt.filter_rows(
        (mt.info.QD >= QD_cutoff) &
        (mt.info.FS < FS_snps_cutoff) &
        (mt.info.MQ > MQ_cutoff) &
        (mt.info.MQRankSum >= MQRankSum_cutoff) &
        (mt.info.ReadPosRankSum >= ReadPosRankSum_cutoff)
    )
    #additional filtering
    mt_snps_filtered2 = mt_snps_filtered.filter_entries(
        (
                mt_snps_filtered.GQ > GQ_cutoff
        )
    )

    return mt_snps_filtered2


def filter_indels_cutoffs(mt:hl.MatrixTable) -> hl.MatrixTable:
    '''
    Filter indels based on Hilary's cutoffs
    :param mt: indels only hail matrixtable
    :return: filtered matrixtable based on the cutoffs listed in cutoffs.py file
    '''
    mt_indels_filtered = mt.filter_rows(

        (mt.info.QD >= QD_cutoff) &
        (mt.info.FS <= FS_indels_cutoff) &
        (mt.info.ReadPosRankSum >= ReadPosRankSum_indels_cutoff)
    )
    return mt_indels_filtered


def join_matrixtables_snps_indels(mt1:hl.MatrixTable, mt2:hl.MatrixTable) -> hl.MatrixTable:
    '''
    Combine the snps and indels matrixtables into one (to be performed after separate filtering for snps and indels)
    :param mt1: the snp only matrix table
    :param mt2: the indel only matrix table
    :return: mt_union: both together
    '''

    mt_union=mt1.union_rows(mt2)
    return mt_union

def overlap_with_file(mt:hl.MatrixTable, bed) -> hl.MatrixTable:
    '''

    :param mt: a matrixtable
    :param bed: the baits file with coordinates for which to filter matrixtable
    :return: a matrixtable with variants overlapping with baits file only
    '''
    baits= hl.import_bed(bed, reference_genome='GRCh38')
    overlapping_mt = mt.filter_rows(hl.is_defined(baits[mt.locus]))
    return overlapping_mt



def calculate_hail_sample_qc(mt:hl.MatrixTable) -> hl.MatrixTable:
    '''
    Compute sample QC metrics
    :param mt: the original matrixtable
    :return: annotate matrixtable with sample qc annotations in column 'sample_qc'
    '''
    mt_with_sampleqc = hl.sample_qc(mt, name='sample_qc')
    return mt_with_sampleqc

def calclulate_hail_variant_qc(mt: hl.MatrixTable) -> hl.MatrixTable:
    '''
    Compute variant qc metrics
    :param mt: the original matrixtable
    :return: annotated matrixtable with variant_qc struct
    '''
    mt_with_variantqc = hl.variant_qc(mt, name='variant_qc')
    return mt_with_variantqc



def missingness_per_individual(mt:hl.MatrixTable) -> hl.MatrixTable:
    '''
    Calclulate missingness per individual and annotate matrixtable with these values
    :param mt: the original matrixtable
    :return: a matrix table annotated with 'missingness' value per sample indexed by column
    '''
    #Call rate per sample
    mt_qc = calculate_hail_sample_qc(mt)
    #missingess = 1 - call_rate
    mt = mt_qc.annotate_cols(missingness_individual = 1 - mt_qc.sample_qc.call_rate )
    return mt


def missingness_per_variant(mt: hl.MatrixTable) -> hl.MatrixTable:
    '''
    Calculate missingess per variant and annotate matrixtable with these values
    :param mt: the original matrixtable
    :return:  the matrixtable with 'missingness' value per variant indexed by row
    '''
    mt_qc=calclulate_hail_variant_qc(mt)
    mt=mt_qc.annotate_rows(missingess_variant = 1 - mt_qc.variant_qc.call_rate)
    return mt

def annotate_samples_with_cohort_info(mt: hl.MatrixTable, cohort_file) -> hl.MatrixTable:
    '''

    :param mt: matrixtable with cohort samples and variants
    :param cohort_file: a txt file with no header line and 2 columns, 1st: for sampleID;  2nd: cohortname; example:/lustre/scratch115/projects/autozyg/new_autozyg_DDD_callset.April2019/sample_list.after_QC.ELGH_BiB_Birm_controls_only.with_cohort_labels.txt
    :return: matrixtable with new column annotation
    '''
    #import the tab delimited file. Note that it is important for joins of tables to have defined keys in the hail tables
    table_cohort = hl.import_table(cohort_file, no_header=True, key='f0')
    #annotate the samples with a new attribute called cohort:
    mt_result = mt.annotate_cols(cohort=table_cohort[mt.s].f1)
    return mt_result


def import_vcf_to_mt(vcf, mt_save_location="None", chr="" ) -> hl.MatrixTable:
    '''

    :param vcf: a bgzipped index vcf file can be choromosome or genomewide
    :param chr: chromosome name if applicable
    :param mt_save_location: s3 location to save the matrixtable for quicker loading in the next steps
    :return: a matrixtable represenation of this vcf
    '''
    if mt_save_location != "None":
        print("Saving to: {}".format(mt_save_location))
        hl.import_vcf(vcf, force_bgz = True).write(mt_save_location+"/chr{}.mt".format(chr), overwrite = True)
        mt=hl.read_matrix_table(mt_save_location+"/chr{}.mt".format(chr))
    else:
        print("Not saving mt to file because no file was provided")
        mt = hl.import_vcf(vcf, force_bgz=True)

    return mt

def import_all_chr_vcf_in_folder(folder, mt_save_location):
    '''

    :param folder: s3 location of all chromosome vcf files
    :param mt_save_location: s3 or local location to save matrix tables per chromosome
    :return:
    '''

    chromosomes=[*range(1,23),"X","Y"]
    #folder='s3://ddd-elgh/vcf'
    for chr in chromosomes:
        vcf="{}/chr{}.vep.vcf.gz".format(folder, chr)
        import_vcf_to_mt(vcf, mt_save_location,chr=chr)

def join_chr_matrixtables_folder(mt_location) -> hl.MatrixTable:
    '''

    :param mt_save_location: folder with chromosome mt files
    :return: combined mt file
    '''
    chr1="1"
    mt_chr1 = hl.read_matrix_table(mt_location+"/chr{}.mt".format(chr1))

    chromosomes = [*range(2, 23), "X", "Y"]
    for chr in chromosomes:
        print("Reading chromosome {} matrix table.".format(chr))
        mt=hl.read_matrix_table(mt_location+"/chr{}.mt".format(chr))
        mt_chr1=mt_chr1.union_rows(mt)


    mt_all= mt_chr1.repartition(n_partitions) #from constants.py - 12800
    mt_all.write(mt_location+"/all_chr.mt", overwrite = True)

    return mt_all

def join_mt_with_vep_json(mt,json)-> hl.MatrixTable:
    '''

    :param mt: matrix table
    :param json: VEP json annotation in gz format
    :return: matrixtable with VEP annotations imported and ready for hail commands
    '''
    vep_annotation = hl.import_table(json, force_bgz= True, no_header=True, types={
        'f0': vep_schema})
    input_str = vep_annotation.f0.input.split("\t")
    vep_annotation = vep_annotation.annotate(in_str=input_str)
    vep_annotation = vep_annotation.key_by(locus=hl.locus(vep_annotation.in_str[0],
                                                          hl.int(vep_annotation.in_str[1])),
                                           alleles=vep_annotation.in_str[3:5])
    mt = mt.annotate_rows(vep=vep_annotation[mt.locus, mt.alleles].f0)
    #mt.entries().select('GT').show()
    #mt.vep.transcript_consequences.lof.describe()
    return mt


def wgs_interval_sample_qc_filter(mt,qc_column = 'sample_qc_and_phenotype', *pass_columns)-> hl.MatrixTable:
    '''

    :param mt: matrixtable with sample QC and annotations
    :param qc_column: name of the column in the matrixtable containing the sample_qc information to filter on
    :param pass_columns: columns to check if they pass filters
    :return: filtered matrixtable after applying sample QC
    Filter on:

    PASS.ID,

    PASS.SampleSwap,

    PASS.Sex,

    PASS.Depth,

    PASS.Median.FreeMix,

    PASS.NRD

    PASS.CHIM

    For association analyses also add

    PASS.DUP=1 and possibly

    PASS.PIHAT=1 (or otherwise we will need to correct for relatedness)


    '''
    mt_filtered = mt.filter_cols(
        (mt[qc_column]['PASS_CHIM'] == 1) & # not in Klaudias list -ask
        (mt[qc_column]['PASS_Depth'] ==1) &
        (mt[qc_column]['PASS_HetHomAlt'] == 1) &
        (mt[qc_column]['PASS_ID'] == 1) &
        (mt[qc_column]['PASS_Median_FreeMix'] == 1) &
        (mt[qc_column]['PASS_NRD'] == 1) &
        (mt[qc_column]['PASS_PIHAT'] == 1) & # for association studies
        (mt[qc_column]['PASS_SampleSwap'] == 1) &
        (mt[qc_column]['PASS_Sex'] == 1) &
        (mt[qc_column]['PASS_DUP'] == 1)  # for association studies
    )
    return mt_filtered



def wgs_interval_variant_qc_filter(mt)-> hl.MatrixTable:
    '''

    :param mt: matrix table with annotations to filter
    :param kwargs: call_rate,vqslod, Hardy Weinberg equillibrium
    :return:
    '''

    mt_filtered = mt.filter_rows(
        (mt.variant_qc.call_rate >= hail_methods.cutoffs.call_rate) &

        (mt.variant_qc.p_value_hwe >= hail_methods.cutoffs.HWE_threshold)
    )

    return mt_filtered

def create_stats(mt) -> hl.Table:
    '''

    :param mt: matrixtable with sample qc and variant qc calculated
    :return: table with only stats for qc to be plotted fot wgs-interval project
    '''
    table_result = mt.select_cols(mt.s, mt.sample_qc.n_snp)

def filter_vqsr(mt, snp_vqsr_table, indel_vqsr_table)-> hl.MatrixTable:
    '''

    :param mt: matrixtable with both snps and indels but having had split multi_alleles
    :param snp_vsqr_table: annotated mt with vqs
    :param indel_vsqr_table:
    :return:
    '''
    #split mt to indels and snp mt
    mt_snps = split_mt_to_snps(mt)
    mt_indels = split_mt_to_indels(mt)

    #annotate each using the snp and indel annotation file respectively
    mt_snps = mt_snps.annotate_rows(vqsr=snp_vqsr_table[mt_snps.locus].VQSLOD)
    mt_indels = mt_indels.annotate_rows(vqsr = indel_vqsr_table[mt_indels.locus].VQSLOD)


    #filter each based on the new column VSQR
    #(mt.vqsr.vqslod >= hail_methods.cutoffs.VQSLOD_score) &
    mt_snps_filtered = mt_snps.filter_rows((mt_snps.vqsr >= SNP_VQSR_threshold))
    mt_indels_filtered = mt_indels.filter_rows((mt_indels.vqsr >=  INDEL_VQSR_threshold))


    #join indel and snp mt together again

    mt_joined = join_matrixtables_snps_indels(mt_snps_filtered, mt_indels_filtered)


    #returned filtered_& joined mt
    return mt_joined




def multi_import_vcf(path) -> hl.MatrixTable:
    '''

    :param path:
    :return:
    '''

    folderfiles = hl.utils.hadoop_ls(path)
    vcf_files= []
    for i in range(0, len(folderfiles)):
        if folderfiles[i]['path'].endswith('all.g.vcf.gz'):
            vcf_files.append(folderfiles[i]['path'])

    #vcf_files=['s3a://interval-wgs-pa10/vcfs/chr1:10001-207666_all.g.vcf.gz',
     #          's3a://interval-wgs-pa10/vcfs/chr1:100231121-100722992_all.g.vcf.gz',
     #          's3a://interval-wgs-pa10/vcfs/chr1:100722993-101221352_all.g.vcf.gz',
      #         's3a://interval-wgs-pa10/vcfs/chr1:101221353-101725726_all.g.vcf.gz',
       #        's3a://interval-wgs-pa10/vcfs/chr1:101725727-102228246_all.g.vcf.gz',
        #       's3a://interval-wgs-pa10/vcfs/chr1:102228247-102451400_all.g.vcf.gz',
         #      's3a://interval-wgs-pa10/vcfs/chr1:10243803-10744548_all.g.vcf.gz',
          #     's3a://interval-wgs-pa10/vcfs/chr1:102451401-102948314_all.g.vcf.gz',
           #    's3a://interval-wgs-pa10/vcfs/chr1:102948315-103446656_all.g.vcf.gz',
            #   's3a://interval-wgs-pa10/vcfs/chr1:103446657-103950770_all.g.vcf.gz'
             #  ]
    #vcf_files = [ f for f in glob.glob(path+'*.vcf.bgz')]
    #pprint(vcf_files)
    mt1 = hl.import_vcf(vcf_files,
                        header_file='s3a://interval-wgs-pa10/vcfs/header_lines.txt',
                        array_elements_required=False,
                        force_bgz=True)

    return mt1
