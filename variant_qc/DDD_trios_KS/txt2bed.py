import csv

un_singl_out="/Users/pa10/Projects/gnomad/trios/forPavlos_100trios_nontransmitted_ac1_syn_variants_2021_07_15_EGA.txt"
trans_singl_out="/Users/pa10/Projects/gnomad/trios/forPavlos_100trios_transmitted_ac2_syn_variants_2021_07_15_EGA.txt"
untr_bed_out="/Users/pa10/Projects/gnomad/trios/forPavlos_100trios_nontransmitted_ac1_b37.bed"
trans_bed_out="/Users/pa10/Projects/gnomad/trios/forPavlos_100trios_transmitted_ac2_b37.bed"
def txt2bed(txt_file, bedfile):
    with open(txt_file, 'r') as csvfile:
        with open(bedfile, 'w') as csvoutput:
            writer = csv.writer(csvoutput,delimiter="\t", lineterminator='\n')
            reader=csv.reader(csvfile, delimiter="\t")
            headers = next(reader)
            all=[]
            headers.insert(2,"stop")
            #print(headers)
            #all.append(headers)
            for row in reader:
                row[0]="chr"+row[0]
                start=int(row[1])-1 
                end=int(row[1])
               
                row[1]=start
                row.insert(2,end)
                #print(row)
                all.append(row[:3])
            
            writer.writerows(all)
            print("Done.")

def main():
    txt2bed(un_singl_out, untr_bed_out)
    txt2bed(trans_singl_out, trans_bed_out)


if __name__=="__main__":
    main()