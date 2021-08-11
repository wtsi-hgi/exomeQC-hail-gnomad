import pandas as pd
import numpy as np
folder = "/Users/pa10/Projects/gnomad"
ancestry_orig = pd.read_csv(f"{folder}/ancestry_data.txt", delimiter="\t")
ancestry_orig.head()
ancestry = ancestry_orig[['eid', '21000-0.0']]
ancestry = ancestry.rename(columns={"21000-0.0": "encoding"})
ancestry.dropna(subset=['encoding'], inplace=True)
ancestry.encoding = ancestry.encoding.astype(int)
encoding_ancestry = pd.read_csv(f"{folder}/coding1001.tsv", delimiter="\t")
df = ancestry.merge(encoding_ancestry, left_on="encoding", right_on="coding")
cohorts = pd.read_csv(f"{folder}/sanger_cohorts.tsv", delimiter="\t")
ukbb_cohorts = cohorts[cohorts['cohort'] == "UKBB"]
sample_ids = df.eid.to_list()
sample_ids2 = ukbb_cohorts.s.to_list()


def find_substrings(df_column, l1):
    # print(df_column)

    for subs in l1:
        if df_column.find(str(subs)) != -1:
            print(f"{df_column}\t{subs}")
            # print(i)
            return(subs)


print("Start running")
ukbb_cohorts['ukbb_sampleid'] = ukbb_cohorts['s'].apply(
    find_substrings, l1=sample_ids)
ukbb_cohorts.to_csv(f"{folder}/matches_samples.tsv", sep="\t", index=False)
