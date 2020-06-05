import re
import csv
map_file = "./UKBB_DDD_ELGH_July19.map"
additional_cohorts = "./EGA_to_study.txt"
samples_list = {}
with open(map_file, newline='') as samples:
    samples_reader = csv.reader(samples, delimiter='\t')
    for s1 in samples_reader:
        if "elgh" in s1[1]:
            samples_list[s1[0]] = "ELGH"
        elif "UK10K" in s1[1] or "ukbb" in s1[1]:
            samples_list[s1[0]] = "UKBB"
        elif "ddd" in s1[1]:
            samples_list[s1[0]] = "DDD"
        elif "chd" in s1[1]:
            samples_list[s1[0]] = "CHD"
        else:
            samples_list[s1[0]] = s1[1]

# print(samples_list)
with open(additional_cohorts, newline='') as samples2:
    samples_reader = csv.reader(samples2, delimiter=' ')
    for s1 in samples_reader:
        samples_list[s1[0]] = s1[1]
        # print(s1[0])

# print(samples_list)
for sample, cohort in samples_list.items():
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
