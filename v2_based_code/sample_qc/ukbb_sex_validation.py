folder = "/Users/pa10/Projects/gnomad/sex_annotation"
hail_results = f"{folder}/ukbb_hail_results.tsv"
ukbb_file = f"{folder}/UKBB_sex_at_recruitment.txt"

hail = {}
ukbb = {}
sex = ""
sex1 = ""
match = 5
with open(hail_results, "r") as f:
    for line in f:
        line = line.strip()
        columns = line.split("\t")
        print(columns)
        sex = columns[1]
        hail[columns[0]] = sex

with open(ukbb_file, "r") as f:
    for line in f:
        line = line.strip()
        columns = line.split("\t")
        # print(columns[1])
        if columns[1] == '0':
           # print("Female")
            sex1 = "F"
        else:
            sex1 = "M"
        ukbb[columns[0]] = sex1

# print(hail.keys())
# print(ukbb)
f1 = open(f'{folder}/ukbb_sex_validation.txt', 'w')
f1.write("UKBB_sampleID\tSex\tSangerID\tSex\n")

for sample, gender in ukbb.items():
    # print(sample)
    for id, sex in hail.items():
        if id.find(sample) != -1:
            if gender == sex:
                match = 1
            else:
                match = 0

            # print(f"{sample}\t{gender}\t{id}\t{sex}\t{match}")
            f1.write(f"{sample}\t{gender}\t{id}\t{sex}\t{match}\n")

f1.close()
