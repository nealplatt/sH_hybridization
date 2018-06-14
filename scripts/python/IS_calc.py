# lift_over_vcf.py <in_bed> <in_vcf> <out_vcf>

import sys

#open necessary files (read in from the cmd line)
vcf_file_in=open(sys.argv[1], 'r')
scores_file_out=open(sys.argv[2], 'w')

vcf_file_in=open(sys.argv[1], 'r')
scores_file_out=open(sys.argv[2], 'w')

#read through the bed file and store the relevant
# information in a dictionary (snp chr and pos)
position = {}
chromosome = {}

for bed_entry in bed_file_in:
    
    snp_id = bed_entry.split("\t")[3].rstrip()

    position[snp_id] = bed_entry.split("\t")[2]
    chromosome[snp_id] = bed_entry.split("\t")[0]

bed_file_in.close()

#alter the chr and start position with the 
#  information from the bed file (stored as a dict)
for vcf_entry in vcf_file_in:

    if vcf_entry.startswith("#"):
        #its a header line
        vcf_file_out.write(vcf_entry)
    else:
        vcf_info = vcf_entry.split("\t")

        snp_id = vcf_entry.split("\t")[2]
    
        if snp_id not in position or snp_id not in chromosome:
            print("FAIL:", snp_id)

        else:
            vcf_info[0] = chromosome[snp_id]
            vcf_info[1] = position[snp_id]
     
            new_vcf_entry = '\t'.join(vcf_info)
            vcf_file_out.write(new_vcf_entry)     


#clean up
vcf_file_in.close()
vcf_file_out.close()
