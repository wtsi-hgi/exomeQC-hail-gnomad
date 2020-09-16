import os
import pandas as pd

'''
A text file with no header line, and one line per sample with the following six fields:

Family ID ('FID')
Within-family ID ('IID'; cannot be '0')
Within-family ID of father ('0' if father isn't in dataset)
Within-family ID of mother ('0' if mother isn't in dataset)
Sex code ('1' = male, '2' = female, '0' = unknown)
Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
'''

trios = "/Users/pa10/Projects/gnomad/trios/trios.txt"
IDmatch = "/Users/pa10/Projects/gnomad/trios/decipher_stable_sanger_ega.txt"

'''
trios
decipher_id	proband_stable_id	mother_stable_id	father_stable_id
294535	DDDP128446	DDDP128444	DDDP128445
259313	DDDP103523	DDDP103522	DDDP103524
IDmatch
decipher_id	person_stable_id	sanger_id	ega_id	is_proband	gender
294535	DDDP128446	DDD_MAIN6090364	EGAN00001351469	1	M
294535:mat	DDDP128444	DDD_MAIN6116947	EGAN00001351470	0	F
294535:pat	DDDP128445	DDD_MAIN6116948	EGAN00001351471	0	M
280693	DDDP119527	DDD_MAIN5749583	EGAN00001301423	1	F
259313	DDDP103523	DDD_MAIN5255949	EGAN00001315854	1	M
# 1. create:
decipher_id	proband_stable_id	mother_stable_id	father_stable_id proband mother father 
294535 DDDP128446	DDDP128444	DDDP128445 EGAN00001351469 EGAN00001351470 EGAN00001351471

2. merge on decipherid 

#merge and create:
#.fam:
# familyid, proband father mother sex phenotype 1 (control)
294535 EGAN00001351469 EGAN00001351471 EGAN00001351470 M 1 
 
'''

df = pd.read_csv(trios, delimiter="\t")
df2 = pd.read_csv(IDmatch, delimiter="\t")
df3 = df.merge(df2, left_on="proband_stable_id", right_on="person_stable_id")
df3 = df3.rename(columns={"ega_id": "proband"})
df3 = df3.drop(columns=["decipher_id_y", "person_stable_id", "sanger_id"])
df3 = df3.rename(columns={"decipher_id_x": "decipher_id"})
df = df3.merge(df2, left_on="mother_stable_id", right_on="person_stable_id")
df = df.rename(columns={"ega_id": "mother"})
df = df.drop(columns=["decipher_id_y", "person_stable_id",
                      "sanger_id", "is_proband_x", "is_proband_y", "gender_y"])
df = df.rename(columns={"decipher_id_x": "decipher_id",
                        "gender_x": "proband_gender"})
df_pat = df.merge(df2, left_on="father_stable_id", right_on="person_stable_id")
df_pat = df_pat.rename(columns={"ega_id": "father"})
df_pat = df_pat.drop(
    columns=["decipher_id_y", "person_stable_id", "sanger_id", "is_proband", "gender"])
df = df_pat.rename(columns={"decipher_id_x": "decipher_id"})
df.to_csv("DDD_trios_sangerIDs.tsv", sep="\t", index=False)
df_fam = df
df_fam["phenotype"] = 1
df_fam = df_fam[["decipher_id", "proband", "father",
                 "mother", "proband_gender", "phenotype"]]
df_fam.to_csv("DDD_trios.fam", sep="\t", header=False, index=False)
