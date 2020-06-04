import re

list1 = ['s3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr10_vcf.1.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr11_vcf.2.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr12_vcf.3.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr13_vcf.4.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr14_vcf.5.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr15_vcf.6.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr16_vcf.7.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr17_vcf.8.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr18_vcf.9.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr19_vcf.10.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr1_vcf.11.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr20_bcftools.vcf.bgz',
         's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr21_vcf.13.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr22_vcf.14.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr2_vcf.15.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr3_vcf.16.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr4_vcf.17.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr5_vcf.18.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr6_vcf.19.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr7_vcf.20.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr8_vcf.21.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chr9_vcf.22.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chrX_vcf.23.vcf.bgz', 's3a://DDD-ELGH-UKBB-exomes/vcfs-chrs/chrY_vcf.24.vcf.bgz']

for vcf in list1:
    m = re.search(r'chr([\d+|\w+]+)_vcf.(\d+).vcf.bgz', vcf)
    if m:
        chromosome = "chr"+m.group(1)
        print(chromosome)
