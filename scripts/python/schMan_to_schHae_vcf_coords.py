# lift_over_vcf.py <in_bed> <in_vcf> <out_vcf>

import sys

#open necessary files (read in from the cmd line)
vcf_file_in=open(sys.argv[1], 'r')
vcf_file_out=open(sys.argv[2], 'w')
#vcf_file_in=open("Smp_127030_cds.vcf", 'r')
#vcf_file_out=open("test.vcf", 'w')

for vcf_entry in vcf_file_in:
    
    if vcf_entry.startswith("#"):
        #its a header line
        vcf_file_out.write(vcf_entry)
    else:
        vcf_info = vcf_entry.split("\t")
        
        #get haem chrom and coord
        chrom = vcf_info[2].split(":")[0]
        coord = vcf_info[2].split(":")[1]
        
        #create new snp_id with mansoni chrom and coord
        id = ':'.join([vcf_info[0], vcf_info[1]])
    
        #get the remainder of the vcf
        remaining='\t'.join(vcf_info[3:])
        new_vcf_entry='\t'.join([chrom, coord, id, remaining])
        
        #now write out
        vcf_file_out.write(new_vcf_entry)

#clean up
vcf_file_in.close()
vcf_file_out.close()
