# build_target_snp_lists.py <ref_panel_bed> <target_list>
import sys

# open necessary files (read in from the cmd line)
bed_file_in = open(sys.argv[1], "r")
snp_list_in = open(sys.argv[2], "r")


target_snps = {}

for snp_entry in snp_list_in:

    snp_id = snp_entry.split("\t")[0].rstrip()

    target_snps[snp_id] = 1

snp_list_in.close()

for bed_entry in bed_file_in:

    snp_id = bed_entry.split("\t")[3].rstrip()
    chromosome = bed_entry.split("\t")[0].rstrip()
    location = bed_entry.split("\t")[2].rstrip()
    target = chromosome + ":" + location

    if snp_id in target_snps:
        print(target)

bed_file_in.close()
