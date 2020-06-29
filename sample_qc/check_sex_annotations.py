import re
import csv
import pandas as pd
# read sex assignment from hail
# make sex M F
# read each of the samples and create a dictionary with sample and sex M F
# end up with
# sample hail recruitment  agreement

# plot if they agree
# list of non agreed
hail_results = "/nfs/users/nfs_m/mercury/pa10/ddd-elgh-ukbb/sex_annotations/chr1_chr20_XY.sex_check.txt.bgz"
DDD = "/nfs/users/nfs_m/mercury/pa10/ddd-elgh-ukbb/sex_annotations/DDD_ega_gender.txt"
# ega_id	gender
# EGAN00001351469	M
# EGAN00001351471	M
# EGAN00001351470	F

elgh = "/nfs/users/nfs_m/mercury/pa10/ddd-elgh-ukbb/sex_annotations/ELGH_samples_with_coverage_and_assigned_sex.txt"
ukbb = "/nfs/users/nfs_m/mercury/pa10/ddd-elgh-ukbb/sex_annotations/UKBB_sex_at_recruitment_EGAN.txt"

df = pd.read_csv(hail_results, compression='gzip', delimiter="\t")
print(df.head())
#                 s is_female    f_stat  n_called  expected_homs  observed_homs
# 0  EGAN00001006259      True  0.140010      1071         743.09            789
# 1  EGAN00001006260     False  0.856640      1071         743.15           1024

df_ddd = pd.read_csv(DDD, delimiter="\t")
print(df_ddd.head())

df_elgh = pd.read_csv(elgh, delimiter="\t")
print(df_elgh.head())

df_ukbb = pd.read_csv(ukbb, delimiter="\t")
print(df_ukbb.head())
