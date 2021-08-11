import hail as hl


def import_vcf_to_mt(vcf, mt_save_location="None", chr="") -> hl.MatrixTable:
    '''

    :param vcf: a bgzipped index vcf file can be choromosome or genomewide
    :param chr: chromosome name if applicable
    :param mt_save_location: s3 location to save the matrixtable for quicker loading in the next steps
    :return: a matrixtable represenation of this vcf
    '''
    if mt_save_location != "None":
        print("Saving to: {}".format(mt_save_location))
        hl.import_vcf(vcf, force_bgz=True).write(
            mt_save_location+"/chr{}.mt".format(chr), overwrite=True)
        mt = hl.read_matrix_table(mt_save_location+"/chr{}.mt".format(chr))
    else:
        print("Not saving mt to file because no file was provided")
        mt = hl.import_vcf(vcf, force_bgz=True)

    return mt


def import_all_chr_vcf_in_folder(folder, mt_save_location):
    '''

    :param folder: s3 location of all chromosome vcf files
    :param mt_save_location: s3 or local location to save matrix tables per chromosome
    :return:
    '''

    chromosomes = [*range(1, 23), "X", "Y"]
    # folder='s3://ddd-elgh/vcf'
    for chr in chromosomes:
        vcf = "{}/chr{}.vep.vcf.gz".format(folder, chr)
        import_vcf_to_mt(vcf, mt_save_location, chr=chr)


def join_chr_matrixtables_folder(mt_location) -> hl.MatrixTable:
    '''

    :param mt_save_location: folder with chromosome mt files
    :return: combined mt file
    '''
    chr1 = "1"
    mt_chr1 = hl.read_matrix_table(mt_location+"/chr{}.mt".format(chr1))

    chromosomes = [*range(2, 23), "X", "Y"]
    for chr in chromosomes:
        print("Reading chromosome {} matrix table.".format(chr))
        mt = hl.read_matrix_table(mt_location+"/chr{}.mt".format(chr))
        mt_chr1 = mt_chr1.union_rows(mt)

    mt_all = mt_chr1.repartition(n_partitions)  # from constants.py - 12800
    mt_all.write(mt_location+"/all_chr.mt", overwrite=True)

    return mt_all
