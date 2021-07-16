import pandas as pd
import csv
import json
from pprint import pprint

trios_names="/Users/pa10/Projects/gnomad/trios/DDD_trios_sangerIDs.tsv"
un_singl="/Users/pa10/Projects/gnomad/trios/forPavlos_100trios_nontransmitted_ac1_syn_variants_2021_07_15.txt"
trans_singl="/Users/pa10/Projects/gnomad/trios/forPavlos_100trios_transmitted_ac2_syn_variants_2021_07_15.txt"
outputfile_json='/Users/pa10/Projects/gnomad/trios/DDD_SangerIDs_matches.json'
un_singl_out="/Users/pa10/Projects/gnomad/trios/forPavlos_100trios_nontransmitted_ac1_syn_variants_2021_07_15_EGA.txt"
trans_singl_out="/Users/pa10/Projects/gnomad/trios/forPavlos_100trios_transmitted_ac2_syn_variants_2021_07_15_EGA.txt"
samples_out="/Users/pa10/Projects/gnomad/trios/forPavlos_100trios_EGA_accessions.txt"
#trios_names=pd.read_csv(trios_names, sep="\t")
#print(trios_names.head)



def create_json_matchIDs(trios_file, outputfile):
    name_matches={}
    with open(trios_names, 'r') as csvfile:
        reader=csv.reader(csvfile,delimiter="\t")
        for row in reader:
        #print(row)
            name_matches[row[1]]=row[4]
            name_matches[row[2]]=row[6]
            name_matches[row[3]]=row[7]

        print(name_matches)
        with open(outputfile, 'w') as outfile:
            json.dump(name_matches, outfile)
        


def match_new_trios(trios_file, json_matches, outputfile):
    dict1={}
    with open(json_matches,'r') as json_file:
        json_dict=json.load(json_file)
        #pprint(json_dict)
    with open(trios_file, 'r') as csvfile:
        with open(outputfile, 'w') as csvoutput:
            writer = csv.writer(csvoutput,delimiter="\t", lineterminator='\n')
            reader=csv.reader(csvfile, delimiter="\t")
            headers = next(reader)[1:]
            all=[]
            headers.append("child_Ega")
            headers.append("mom_Ega")
            headers.append("dad_Ega")
            all.append(headers)
            for row in reader:
            #dict1[row[0]] = {key: value for key, value in zip(headers, row[1:])}    
                if row[20] in json_dict:
                    child_ega=json_dict[row[20]]
                    row.append(child_ega)
               
                if row[21] in json_dict:
                    mom_ega=json_dict[row[21]]
                    row.append(mom_ega)
                
                if row[22] in json_dict:
                    dad_ega=json_dict[row[22]]
                    row.append(dad_ega)
                else:
                    print("Not in dictionary")

                all.append(row)
            #19,20,21,22

            writer.writerows(all)
                #print(row[21])
                #print(row[22])
       # pprint(dict1)

def print_EGA_accessions(singletons_file):
    accessions=[]
    with open(singletons_file, 'r') as cvsfile:
            
            reader=csv.reader(cvsfile, delimiter="\t")
            headers = next(reader)[1:]
            for row in reader:
                accessions.append(row[28])
                accessions.append(row[29])
                accessions.append(row[30])
    
    return accessions
               

def main():
    # Step 1 done
    #create_json_matchIDs(trios_names,outputfile_json)

    #match_new_trios(un_singl, outputfile_json, un_singl_out)
    #match_new_trios(trans_singl,outputfile_json,trans_singl_out)

   accessions_list1= print_EGA_accessions(un_singl_out)
   accessions_list2=print_EGA_accessions(trans_singl_out)
   accessions_list=accessions_list1 + accessions_list2
   accessions_list=sorted(set(accessions_list))
   with open(samples_out, 'w') as outfile:
       outfile.writelines(["%s\n" % item for item in accessions_list])

if __name__=="__main__":
    main()