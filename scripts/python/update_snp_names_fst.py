# update_snp_names_fst.py <fst_table> <in_vcf>

import sys
import re

# open necessary files (read in from the cmd line)
fst_in = open(sys.argv[1], "r")
vcf_in = open(sys.argv[2], "r")


# read through the bed file and store the relevant
# information in a dictionary (snp chr and pos)
vcf_snp_ids = {}

for vcf_entry in vcf_in:

    if not vcf_entry.startswith("#"):
        chrom = vcf_entry.split("\t")[0]
        pos = vcf_entry.split("\t")[1]
        fst_id = chrom + "\t" + pos

        vcf_snp_ids[fst_id] = vcf_entry.split("\t")[2]
vcf_in.close()

for fst_entry in fst_in:
    fst_entry = fst_entry.rstrip()

    if not fst_entry.startswith("CHROM"):
        fst_id = fst_entry.split("\t")[0] + "\t" + fst_entry.split("\t")[1]
        # vcf_snp_id = vcf_snp_ids[fst_id]

        print(vcf_snp_ids[fst_id] + "\t" + fst_entry)
        # print(vcf_snp_id)
    else:
        print("SNP_NAME" + "\t" + fst_entry)
fst_in.close()
