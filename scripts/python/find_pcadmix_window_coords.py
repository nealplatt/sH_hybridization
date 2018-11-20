import sys

#open necessary files (read in from the cmd line)
#windows_in=open(sys.argv[1], 'r')
#vcf_file_in=open(sys.argv[2], 'r')


windows_in=open("pcadmix_window.start_stop", 'r')
vcf_file_in=open("../build_snp_panel/cohort_snps_schMan.vcf", 'r')


pos={}
chr={}

for vcf_entry in vcf_file_in:
    
    if not vcf_entry.startswith("#"):
        vcf_info = vcf_entry.split("\t")
        
        snp_id = vcf_info[2]
        snp_pos = vcf_info[1]
        snp_chr = vcf_info[0]
        
        pos[snp_id]=snp_pos
        chr[snp_id]=snp_chr

for window_entry in windows_in: 
    marker_info = window_entry.split()
    window_name = marker_info[0]
    start_snp = marker_info[1]
    stop_snp = marker_info[-1]
        
    window_start_chr=chr[start_snp]
    window_start_pos=pos[start_snp]
    window_end_chr=chr[stop_snp]
    window_end_pos=pos[stop_snp]

    print(window_name, window_start_chr, window_start_pos, window_end_chr_window_end_pos)


#clean up
windows_in.close()
vcf_file_in.close()

