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

df = pd.read_csv(hail_results, compression='gzip', delimiter="\t")
print(df.head())
