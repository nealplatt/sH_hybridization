import sys

#open necessary files (read in from the cmd line)
snps_in=open(sys.argv[1], 'r')
windows_coords_out=open(sys.argv[2], 'r')


pos={}
chr={}
snps_in=open("niger_tz-out_auto_sprime_2018-08-13.score", 'r')
windows_coords_out=open("niger_tz-out_auto_sprime_2018-08-13.cul_pos", 'w')

SM_V7_1_len=88881357
SM_V7_2_len=48130368
SM_V7_3_len=50458499
SM_V7_4_len=47279781
SM_V7_5_len=25256119
SM_V7_6_len=24989083
SM_V7_7_len=19288021


chr_length={}
chr_length["SM_V7_1"]=0
chr_length["SM_V7_2"]=chr_length["SM_V7_1"] + SM_V7_1_len
chr_length["SM_V7_3"]=chr_length["SM_V7_2"] + SM_V7_2_len
chr_length["SM_V7_4"]=chr_length["SM_V7_3"] + SM_V7_3_len
chr_length["SM_V7_5"]=chr_length["SM_V7_4"] + SM_V7_4_len
chr_length["SM_V7_6"]=chr_length["SM_V7_5"] + SM_V7_5_len
chr_length["SM_V7_7"]=chr_length["SM_V7_6"] + SM_V7_6_len

i=0

for snp_entry in snps_in: 
    
    if i >0:
        marker_info = snp_entry.split()
        chr = marker_info[0]
        pos = int(marker_info[1])
        score = marker_info[7]
        
        cul_pos=(pos)+chr_length[chr]
        
        entry='\t'.join([str(chr), str(pos), str(cul_pos), str(score)])
        
        windows_coords_out.write(entry + '\n')    
    
    i=1
    

#clean up
snps_in.close()
windows_coords_out.close()
