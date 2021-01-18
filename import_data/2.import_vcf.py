
# Pavlos Antoniou
# pa10@sanger.ac.uk
# 18/01/2020
# Import vcf files from s3 location to hail and convert to matrixtable

import os
import hail as hl
import pyspark
import json
import sys
import re
from pathlib import Path

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

partitions = 1000

# Create a file with only the header of the VCFs. Helps hail process import faster.
# Define location of input VCFs and their header file
s3location_input = storage["ddd_elgh_ukbb_exomes"]["s3"]["vcfs"]
vcf_header_file = storage["ddd_elgh_ukbb_exomes"]["s3"]["header_file"]


if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory
    tmp_dir = "hdfs://spark-master:9820/"
    temp_dir = os.path.join(os.environ["HAIL_HOME"], "tmp")
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])

    #####################################################################
    ###################### INPUT DATA  ##############################
    #####################################################################
    objects = hl.utils.hadoop_ls(s3location_input)
    print("Reading vcf files")
    vcfs = [vcf["path"] for vcf in objects if vcf["path"].endswith(".bgz")]
    print(vcfs)
    for vcf in vcfs:
        print(vcf)
        m = re.search(r'chr([\d+|\w+]+)_vcf.(\d+).vcf.bgz', vcf)

        if m:
            chromosome = "chr"+m.group(1)
            print(chromosome)
            mt = hl.import_vcf(vcf, array_elements_required=False,
                               force_bgz=True, header_file=vcf_header_file)

            print("Imported vcf file for" + chromosome)

            print("Write to disk:")

            mt.write(
                f"{tmp_dir}/ddd-elgh-ukbb/{chromosome}.mt", overwrite=True)
            print(f"Wrote matrixtable for {chromosome}")
