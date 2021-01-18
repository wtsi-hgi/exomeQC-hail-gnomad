#!/usr/bin/env bash


 # while read -r FILE; do COUNT=$(( $COUNT + 1 )); gs_name="$(basename "${FILE}" ".gz").${COUNT}.bgz"; s3cmd sync --multipart-chunk-size-mb=100 "${FILE}" "${BUCKET}/${gs_name}"; echo ${COUNT}; done < "${CHUNK}"

declare CHUNK="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/output_vcf/s3_gnomad_upload.fofn"
declare BUCKET="s3://DDD-ELGH-UKBB-exomes/vcfs"
while read -r FILE; do
	COUNT=$(( $COUNT + 1 ))
	gs_name="$(basename "${FILE}" ".gz").${COUNT}.bgz"
	#s3cmd put "${FILE}" "${BUCKET}/${CHROMOSOME}-vcf/${gs_name}"
	s3cmd sync --multipart-chunk-size-mb=100 "${FILE}" "${BUCKET}/${gs_name}"
	echo ${COUNT}
done < "${CHUNK}"



declare CHUNK="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/output_vcf/stripped_vcf/vep/veps.fofn"
declare BUCKET="s3://DDD-ELGH-UKBB-exomes/vep"
while read -r FILE; do
	COUNT=$(( $COUNT + 1 ))
	gs_name="$(basename "${FILE}" ".gz").bgz"
	#s3cmd put "${FILE}" "${BUCKET}/${CHROMOSOME}-vcf/${gs_name}"
	s3cmd sync --multipart-chunk-size-mb=100 "${FILE}" "${BUCKET}/${gs_name}"
	echo ${COUNT}
done < "${CHUNK}"
#########

cd /software/common-apps/bcftools-1.9-220/bin
declare CHUNK="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/output_vcf/stripped_vcf/vep/veps.fofn"
declare output="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/output_vcf/stripped_vcf/vep/vep_consequences_worst_synonymous.tsv"
while read -r FILE; do
echo $FILE
bcftools +split-vep $FILE -f '%CHROM\t%POS\t%REF\t%ALT\t%Consequence\n' -s worst:synonymous >> "$output"

done < "${CHUNK}"

############################
declare CHUNK="/lustre/scratch118/humgen/resources/gnomAD/release-3.1.1/gnomad_sites.fofn"
declare BUCKET="s3://gnomad-release-3.1.1"
while read -r FILE; do
echo $FILE
gs_name="$(basename "${FILE}")"
echo $gs_name
s3cmd sync --multipart-chunk-size-mb=100 "${FILE}" "${BUCKET}/${gs_name}"
done < "${CHUNK}"

declare CHUNK="/lustre/scratch118/humgen/resources/gnomAD/release-3.1.1/gnomad_sites.fofn"
declare output="/lustre/scratch118/humgen/resources/gnomAD/release-3.1.1/sites_AF.tsv"
while read -r FILE; do
echo $FILE
bcftools query  -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' $FILE >> "$output"
done < "${CHUNK}"


#bcftools query  -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' /lustre/scratch118/humgen/resources/gnomAD/release-3.0/gnomad.genomes.r3.0.sites.vcf.bgz -o /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gnomad_3.0_sites_AF.tsv