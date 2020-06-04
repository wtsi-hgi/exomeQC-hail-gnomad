import re
import csv
map_file = "./UKBB_DDD_ELGH_July19.map"


with open(map_file, newline='') as samples:
    samples_reader = csv.reader(samples, delimiter='\t')
    for s1 in samples_reader:
        if "elgh" in s1[1]:
            print(f'{s1[0]}\tELGH')
        elif "UK10K" in s1[1] or "ukbb" in s1[1]:
            print(f'{s1[0]}\tUKBB')
        elif "ddd" in s1[1]:
            print(f'{s1[0]}\tDDD')
        elif "chd" in s1[1]:
            print(f'{s1[0]}\tCHD')
        else:
            print(f"UNKOWN_STUDY\t{s1[0]}\t{s1[1]}")
