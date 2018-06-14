# lift_over_vcf.py <in_bed> <in_vcf> <out_vcf>

import sys

#open necessary files (read in from the cmd line)
markers_in=open(sys.argv[1], 'r')
vcf_file_in=open(sys.argv[2], 'r')
windows_coords_out=open(sys.argv[2], 'r')

pos={}
chr={}
markers_in=open("bovis_vs_tz.markers.start_end", 'r')
vcf_file_in=open("../beagle_impute/cohort_snps_schMan_autosomal_beagle.vcf", 'r')
windows_coords_out=open("test", 'w')
for vcf_entry in vcf_file_in:
    
    if not vcf_entry.startswith("#"):
        vcf_info = vcf_entry.split("\t")
        
        snp_id = vcf_info[2]
        snp_pos = vcf_info[1]
        snp_chr = vcf_info[0]
        
        pos[snp_id]=snp_pos
        chr[snp_id]=snp_chr

for window_entry in markers_in: 
    marker_info = window_entry.split()
    window_name = marker_info[0]
    start_snp = marker_info[1]
    stop_snp = marker_info[2]
    
    entry='\t'.join([window_name, chr[start_snp], pos[start_snp], chr[stop_snp], pos[stop_snp]])
    
    windows_coords_out.write(entry + '\n')    


#clean up
markers_in.close()
vcf_file_in.close()
windows_coords_out.close()
