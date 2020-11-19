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