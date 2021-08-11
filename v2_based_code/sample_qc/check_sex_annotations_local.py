import re
import csv
import pandas as pd
folder = "/Users/pa10/Projects/gnomad/sex_annotation"
hail_results = f"{folder}/chr1_chr20_XY.sex_check.txt.bgz"
cohorts = f"{folder}/new_cohorts.tsv"
DDD = f"{folder}/DDD_ega_gender.txt"
elgh = f"{folder}/ELGH_samples_with_coverage_and_assigned_sex.txt"
ukbb = f"{folder}/UKBB_sex_at_recruitment.txt"


df = pd.read_csv(hail_results, compression='gzip', delimiter="\t")
print(df.head())

# cohort info
cohorts_df = pd.read_csv(cohorts, delimiter="\t")
cohorts_df.head()

df1 = pd.merge(df, cohorts_df, how='left', left_on="s",
               right_on="sample", left_index=True)

df1.head()

df = df1[["s", "is_female", "cohort"]]

df.head()

df_ddd = pd.read_csv(DDD, delimiter="\t")
df_ddd.head()

df2 = pd.merge(df, df_ddd, how='left', left_on="s",
               right_on="ega_id", left_index=True)

df2.head()

df2 = df2.rename(columns={"gender": "DDD_gender"})
df2 = df2.drop(columns=['ega_id'])

DDDdf = df2[df2["cohort"] == "DDD"]
DDDdf.head()

df2.head()

df_elgh = pd.read_csv(elgh, delimiter="\t")
df_elgh.head()

df_elgh = df_elgh.drop(columns=["sample", "bam", "total", "chrX", "chrY",
                                "batch", "SupplierID", "Sex", "blah", "low.coverage", "Y.to.X"])
df_elgh = df_elgh.rename(columns={'sex.assigned': 'ELGH_gender'})
df_elgh.loc[df_elgh['ELGH_gender'] == "female", ['ELGH_gender']] = 'F'
df_elgh.loc[df_elgh['ELGH_gender'] == "male", ['ELGH_gender']] = 'M'
df3 = pd.merge(df2, df_elgh, how='left', left_on="s",
               right_on="EGAID", left_index=True)

df_elgh.head()

df3.head()

df3.drop(columns="EGAID")

cohorts_list = df3.cohort.unique()
cohorts_list


df3.loc[df3['is_female'] == True, ['is_female']] = 'F'
df3.loc[df3['is_female'] == False, ['is_female']] = 'M'
df3 = df3.rename(columns={'is_female': 'hail_result'})
df3.head()

df3.head()


df3.to_csv(f"{folder}/validation_sex_check_DDD_ELGH.tsv",
           sep="\t", index=False)

ukbbdf = df3[df3["cohort"] == "UKBB"]
ukbbdf

df_ukbb = pd.read_csv(ukbb, delimiter="\t")
df_ukbb
df_ukbb = df_ukbb.rename(columns={'eid': 'sample', '31-0.0': 'ukbb_gender'})

df_ukbb.loc[df_ukbb['ukbb_gender'] == 1, ['ukbb_gender']] = 'M'
df_ukbb.loc[df_ukbb['ukbb_gender'] == 0, ['ukbb_gender']] = 'F'


df_ukbb.head()

samples_list = ukbbdf['s'].tolist()
minisample_list = df_ukbb['sample'].tolist()
results = []
results = []
ukb_indices = df_ukbb
j = 0
for sample in minisample_list:
    res = [i for i in samples_list if str(sample) in i]
    results.append(sample)

    j += 1
    print(j)
with open(f'{folder}/ukbb_sex_validation_samples_overlap.text', 'w') as f:
    for item in results:
        f.write(str(item))

# per cohort
# new column match True false if same sex
# plot matched non-matched per cohort and percentage
# save figure
