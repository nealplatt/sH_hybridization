# vcf_to_map_constant_recomb.py <BP = 1CM> <in_vcf>

import sys

# open necessary files (read in from the cmd line)
recomb_rate = float(sys.argv[1])
# vcf_file_in = open(sys.argv[2], 'r')
vcf_file_in = open("clean_auto.vcf", "r")

chromosomes = []
positions = []

for vcf_entry in vcf_file_in:

    if not vcf_entry.startswith("#"):
        chromosomes.append(vcf_entry.split("\t")[0])
        positions.append(vcf_entry.split("\t")[1])

vcf_file_in.close()


for i in range(0, len(chromosomes)):

    # cycle through snps
    if i + 1 < len(chromosomes):
        current_chr = chromosomes[i]
        current_pos = int(positions[i])
        next_chr = chromosomes[i + 1]
        next_pos = int(positions[i + 1])

        # set the start pos on the first pas through list of snps
        if i == 0:
            start_pos = current_pos

        # if snps are on teh same chromosome
        if current_chr == next_chr:
            # calculate recomb rate form the first snp on the chromosome
            cm_distance = "%f" % ((current_pos - start_pos) / recomb_rate)
        # if the snps are on different chromosomes
        else:
            # the recomb rate between the two is NA
            cm_distance = "%f" % ((current_pos - start_pos) / recomb_rate)
            # and the start pos for subsequent snps is calculated from the 1st
            #  pos on the next chromosome
            start_pos = int(positions[i + 1])

        print(current_chr, "\t.\t", cm_distance, "\t", current_pos)
