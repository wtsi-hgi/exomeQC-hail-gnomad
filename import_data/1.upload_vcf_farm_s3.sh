#!/usr/bin/env bash


 # while read -r FILE; do COUNT=$(( $COUNT + 1 )); gs_name="$(basename "${FILE}" ".gz").${COUNT}.bgz"; s3cmd sync --multipart-chunk-size-mb=100 "${FILE}" "${BUCKET}/${gs_name}"; echo ${COUNT}; done < "${CHUNK}"

declare CHUNK="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/output_vcf/ddd-elgh-ukbb-vcfs.fofn"
declare BUCKET="s3://DDD-ELGH-UKBB-exomes/vcfs"
while read -r FILE; do
	COUNT=$(( $COUNT + 1 ))
	gs_name="$(basename "${FILE}" ".gz").${COUNT}.bgz"
	#s3cmd put "${FILE}" "${BUCKET}/${CHROMOSOME}-vcf/${gs_name}"
	s3cmd sync --multipart-chunk-size-mb=100 "${FILE}" "${BUCKET}/${gs_name}"
	echo ${COUNT}
done < "${CHUNK}"
