#!/usr/bin/env bash
# Pavlos Antoniou
# pa10@sanger.ac.uk
# 18/01/2020

# Copy files from farm to s3 using s3cmd
# Location of fofn file with the paths to each vcf file -change to your path
declare CHUNK="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/output_vcf/s3_gnomad_upload.fofn"
# change to path of the bucket on s3 
declare BUCKET="s3://DDD-ELGH-UKBB-exomes/vcfs"
while read -r FILE; do
	COUNT=$(( $COUNT + 1 ))
	gs_name="$(basename "${FILE}" ".gz").${COUNT}.bgz"
	s3cmd sync --multipart-chunk-size-mb=100 "${FILE}" "${BUCKET}/${gs_name}"
	echo ${COUNT}
done < "${CHUNK}"

