import re
import csv
map_file = "/nfs/users/nfs_m/mercury/pa10/ddd-elgh-ukbb/UKBB_DDD_ELGH_July19.map"

map_file2 = "/nfs/users/nfs_m/mercury/pa10/ddd-elgh-ukbb/sample_list.before_QC.with_cohort_labels.txt"

map_file3 = "/nfs/users/nfs_m/mercury/pa10/ddd-elgh-ukbb/IHTP_ISC_British_Autozygosity_Populations_Resource_Part_5.all_IDs.180119.txt"
# 15001608192016	sc_autozygELGH7220216	EGAN00001833179	Female	NULL
map_file4 = "/nfs/users/nfs_m/mercury/pa10/ddd-elgh-ukbb/EGA_to_study.txt"

# First read the samples
samples_list = {}
ELGH_1 = "ELGH_Parts1-4"
ELGH_5 = "ELGH_Part5"

with open(map_file, newline='') as samples:
    samples_reader = csv.reader(samples, delimiter='\t')
    for s1 in samples_reader:
        if "elgh" in s1[1]:
            samples_list[s1[0]] = ELGH_1
        elif "UK10K" in s1[1] or "ukbb" in s1[1]:
            samples_list[s1[0]] = "UKBB"
        elif "ddd" in s1[1]:
            samples_list[s1[0]] = "DDD"
        elif "chd" in s1[1]:
            samples_list[s1[0]] = "CHD"
        else:
            samples_list[s1[0]] = s1[1]

# print(samples_list)
with open(map_file2, newline='') as samples2:
    samples_reader = csv.reader(samples2, delimiter='\t')

    for s1 in samples_reader:
        if s1[0] in samples_list.keys():
            samples_list[s1[0]] = s1[1]
        # print(s1[0])
with open(map_file3, newline='') as samples3:
    samples_reader = csv.reader(samples3, delimiter='\t')
    for s1 in samples_reader:
        if s1[2] in samples_list.keys():
            samples_list[s1[2]] = ELGH_5

with open(map_file4, newline='') as samples4:
    samples_reader = csv.reader(samples4, delimiter=' ')
    for s1 in samples_reader:
        if s1[0] in samples_list.keys() and "ELGH" not in samples_list[s1[0]]:
            if "IHTP_ISC_British" in s1[2] or "SEQCAP_WES_AUTOZYGOSITY" in s1[2]:
                samples_list[s1[0]] = ELGH_1
            else:
                samples_list[s1[0]] = s1[2]

# print(samples_list)
print(f"sample\tcohort")
for sample, cohort in samples_list.items():
    if cohort == "ELGH":
        cohort = ELGH_1
    print(f"{sample}\t{cohort}")
# if "elgh" in s1[1]:
#     print(f'{s1[0]}\tELGH')
# elif "UK10K" in s1[1] or "ukbb" in s1[1]:
#     print(f'{s1[0]}\tUKBB')
# elif "ddd" in s1[1]:
#     print(f'{s1[0]}\tDDD')
# elif "chd" in s1[1]:
#     print(f'{s1[0]}\tCHD')
# else:
#     print(f"UNKOWN_STUDY\t{s1[0]}\t{s1[1]}")