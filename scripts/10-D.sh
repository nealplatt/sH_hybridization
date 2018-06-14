#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir D

cd D


#in python

import allel
import numpy as np  
import sys
import itertools

#read in the vcf data
callset=allel.read_vcf('../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf', log=sys.stdout)
gt=allel.GenotypeArray(callset['calldata/GT'])

#get allele counts for each locus
ac=gt.count_alleles()


#designate the population (index) in the allele count array
egypt_pop=[0,1,9] 
bov_pop=2
mat_pop=[3,7,8]
inter_pop=4
curs_pop=5
marg_pop=6      
niger_pop=list(range(10,58))                                                            
tz_pop=list(range(58,104))

bov_group_pop=[bov_pop, curs_pop, marg_pop]
haem_group_pop=egypt_pop + niger_pop + tz_pop

#now count the alleles in the array for each pop
egypt_ac=gt.count_alleles(subpop=egypt_pop)
bov_ac=gt.count_alleles(subpop=[bov_pop])
mat_ac=gt.count_alleles(subpop=mat_pop)
inter_ac=gt.count_alleles(subpop=[inter_pop])
curs_ac=gt.count_alleles(subpop=[curs_pop])
marg_ac=gt.count_alleles(subpop=[marg_pop])
niger_ac=gt.count_alleles(subpop=niger_pop)
tz_ac=gt.count_alleles(subpop=tz_pop)

bov_group_ac=gt.count_alleles(subpop=bov_group_pop)
haem_group_ac=gt.count_alleles(subpop=haem_group_pop)

#calculate D (genome wide)

allel.average_patterson_d(niger_ac, tz_ac, bov_ac, marg_ac, blen=100)
#D=0.5446731600901533
#SE=0.020470047021376385
#Z=26.60830038745705

allel.average_patterson_d(niger_ac, tz_ac, curs_ac, marg_ac, blen=100)
#D=0.39719782215941596
#SE=0.02204818857784749
#Z=18.014986617018195

#D sliding window (niger vs. tz)----------------------------------------------
#identify where in the array each chr starts and begins
#...starts
chr1_end_idx=np.amax(np.where(callset['variants/CHROM']=="SM_V7_1"))
chr2_end_idx=np.amax(np.where(callset['variants/CHROM']=="SM_V7_2"))
chr3_end_idx=np.amax(np.where(callset['variants/CHROM']=="SM_V7_3"))
chr4_end_idx=np.amax(np.where(callset['variants/CHROM']=="SM_V7_4"))
chr5_end_idx=np.amax(np.where(callset['variants/CHROM']=="SM_V7_5"))
chr6_end_idx=np.amax(np.where(callset['variants/CHROM']=="SM_V7_6"))
chr7_end_idx=np.amax(np.where(callset['variants/CHROM']=="SM_V7_7"))
#...ends
chr1_start_idx=np.amin(np.where(callset['variants/CHROM']=="SM_V7_1"))
chr2_start_idx=np.amin(np.where(callset['variants/CHROM']=="SM_V7_2"))
chr3_start_idx=np.amin(np.where(callset['variants/CHROM']=="SM_V7_3"))
chr4_start_idx=np.amin(np.where(callset['variants/CHROM']=="SM_V7_4"))
chr5_start_idx=np.amin(np.where(callset['variants/CHROM']=="SM_V7_5"))
chr6_start_idx=np.amin(np.where(callset['variants/CHROM']=="SM_V7_6"))
chr7_start_idx=np.amin(np.where(callset['variants/CHROM']=="SM_V7_7"))

window_len=500
step=500

chr1_d=allel.moving_patterson_d(niger_ac, tz_ac, bov_ac, marg_ac, size=window_len, step=step, start=chr1_start_idx, stop=chr1_end_idx) 
chr2_d=allel.moving_patterson_d(niger_ac, tz_ac, bov_ac, marg_ac, size=window_len, step=step, start=chr2_start_idx, stop=chr2_end_idx) 
chr3_d=allel.moving_patterson_d(niger_ac, tz_ac, bov_ac, marg_ac, size=window_len, step=step, start=chr3_start_idx, stop=chr3_end_idx) 
chr4_d=allel.moving_patterson_d(niger_ac, tz_ac, bov_ac, marg_ac, size=window_len, step=step, start=chr4_start_idx, stop=chr4_end_idx) 
chr5_d=allel.moving_patterson_d(niger_ac, tz_ac, bov_ac, marg_ac, size=window_len, step=step, start=chr5_start_idx, stop=chr5_end_idx) 
chr6_d=allel.moving_patterson_d(niger_ac, tz_ac, bov_ac, marg_ac, size=window_len, step=step, start=chr6_start_idx, stop=chr6_end_idx) 
chr7_d=allel.moving_patterson_d(niger_ac, tz_ac, bov_ac, marg_ac, size=window_len, step=step, start=chr7_start_idx, stop=chr7_end_idx) 

#chr averages
np.average(chr1_d[~np.isnan(chr1_d)])
np.average(chr2_d[~np.isnan(chr2_d)]) 
np.average(chr3_d[~np.isnan(chr3_d)]) 
np.average(chr4_d[~np.isnan(chr4_d)]) 
np.average(chr5_d[~np.isnan(chr5_d)]) 
np.average(chr6_d[~np.isnan(chr6_d)]) 
np.average(chr7_d[~np.isnan(chr7_d)]) 

#chr1 0.396
#chr2 0.404
#chr3 0.385
#chr4 0.559
#chr5 0.667
#chr6 0.598
#chr7 0.312

len(chr1_d) #714
len(chr2_d) #371
len(chr3_d) #374
len(chr4_d) #364
len(chr5_d) #147
len(chr6_d) #168
len(chr7_d) #112
#2250

pos=list(callset['variants/POS'])

chr1_start=pos[chr1_start_idx: (chr1_end_idx-window_len): step]
chr2_start=pos[chr2_start_idx: (chr2_end_idx-window_len): step]
chr3_start=pos[chr3_start_idx: (chr3_end_idx-window_len): step]
chr4_start=pos[chr4_start_idx: (chr4_end_idx-window_len): step]
chr5_start=pos[chr5_start_idx: (chr5_end_idx-window_len): step]
chr6_start=pos[chr6_start_idx: (chr6_end_idx-window_len): step]
chr7_start=pos[chr7_start_idx: (chr7_end_idx-window_len): step]

chr1_end=pos[(chr1_start_idx+window_len): chr1_end_idx: step]
chr2_end=pos[(chr2_start_idx+window_len): chr2_end_idx: step]
chr3_end=pos[(chr3_start_idx+window_len): chr3_end_idx: step]
chr4_end=pos[(chr4_start_idx+window_len): chr4_end_idx: step]
chr5_end=pos[(chr5_start_idx+window_len): chr5_end_idx: step]
chr6_end=pos[(chr6_start_idx+window_len): chr6_end_idx: step]
chr7_end=pos[(chr7_start_idx+window_len): chr7_end_idx: step]

chr1_window_pos=list(np.average(np.array((chr1_start, chr1_end)), axis=0))
chr2_window_pos=list(np.average(np.array((chr2_start, chr2_end)), axis=0))
chr3_window_pos=list(np.average(np.array((chr3_start, chr3_end)), axis=0))
chr4_window_pos=list(np.average(np.array((chr4_start, chr4_end)), axis=0))
chr5_window_pos=list(np.average(np.array((chr5_start, chr5_end)), axis=0))
chr6_window_pos=list(np.average(np.array((chr6_start, chr6_end)), axis=0))
chr7_window_pos=list(np.average(np.array((chr7_start, chr7_end)), axis=0))

#numbers from the VCF ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf
chr1_len=88881357
chr2_len=48130368
chr3_len=50458499
chr4_len=47279781
chr5_len=25256119
chr6_len=24989083
chr7_len=19288021

start=0
chr1_cul_pos=[x+0 for x in chr1_window_pos]
start=chr1_len
chr2_cul_pos=[x+start for x in chr2_window_pos]
start+=chr2_len
chr3_cul_pos=[x+start for x in chr3_window_pos]
start+=chr3_len
chr4_cul_pos=[x+start for x in chr4_window_pos]
start+=chr4_len
chr5_cul_pos=[x+start for x in chr5_window_pos]
start+=chr5_len
chr6_cul_pos=[x+start for x in chr6_window_pos]
start+=chr6_len
chr7_cul_pos=[x+start for x in chr7_window_pos]

cul_pos=""
cul_pos=chr1_cul_pos
cul_pos+=chr2_cul_pos
cul_pos+=chr3_cul_pos
cul_pos+=chr4_cul_pos
cul_pos+=chr5_cul_pos
cul_pos+=chr6_cul_pos
cul_pos+=chr7_cul_pos

chr=""
chr=list(itertools.repeat("SM_V7_1", len(chr1_cul_pos)))
chr+=list(itertools.repeat("SM_V7_2", len(chr2_cul_pos)))
chr+=list(itertools.repeat("SM_V7_3", len(chr3_cul_pos)))
chr+=list(itertools.repeat("SM_V7_4", len(chr4_cul_pos)))
chr+=list(itertools.repeat("SM_V7_5", len(chr5_cul_pos)))
chr+=list(itertools.repeat("SM_V7_6", len(chr6_cul_pos)))
chr+=list(itertools.repeat("SM_V7_7", len(chr7_cul_pos)))

chr=""
chr=list(itertools.repeat("SM_V7_1", 214))
chr+=list(itertools.repeat("SM_V7_2", 111))
chr+=list(itertools.repeat("SM_V7_3", 112))
chr+=list(itertools.repeat("SM_V7_4", 109))
chr+=list(itertools.repeat("SM_V7_5", 44))
chr+=list(itertools.repeat("SM_V7_6", 50))
chr+=list(itertools.repeat("SM_V7_7", 34))


D_per_window=list(chr1_d) + list(chr2_d) + list(chr3_d) + list(chr4_d) + list(chr5_d) + list(chr6_d) + list(chr7_d) 
window_pos=chr1_window_pos + chr2_window_pos + chr3_window_pos + chr4_window_pos + chr5_window_pos + chr6_window_pos + chr7_window_pos


window_d=np.column_stack((chr, cul_pos, window_pos, D_per_window))

# remove if count is lt 10 #######

#save to file
np.savetxt("window_D_500snp-500snp_autosomal.csv", window_d, fmt='%s', delimiter=",")

#now plot this in R so that each chromosome has a different color



##13536
#chr1_cul_pos=[x+0 for x in chr1_window_pos]
#start=chr1_len
#chr2_cul_pos=[x+start for x in chr2_window_pos]
#start+=chr2_len
#chr3_cul_pos=[x+start for x in chr3_window_pos]
#start+=chr3_len
#chr4_cul_pos=[x+start for x in chr4_window_pos]
#start+=chr4_len
#chr5_cul_pos=[x+start for x in chr5_window_pos]
#start+=chr5_len
#chr6_cul_pos=[x+start for x in chr6_window_pos]
#start+=chr6_len
#chr7_cul_pos=[x+start for x in chr7_window_pos]


#cul_pos=chr1_cul_pos
#cul_pos+=chr2_cul_pos
#cul_pos+=chr3_cul_pos
#cul_pos+=chr4_cul_pos
#cul_pos+=chr5_cul_pos
#cul_pos+=chr6_cul_pos
#cul_pos+=chr7_cul_pos

#chr=""
#chr=list(itertools.repeat("SM_V7_1", len(chr1_cul_pos)))
#chr+=list(itertools.repeat("SM_V7_2", len(chr2_cul_pos)))
#chr+=list(itertools.repeat("SM_V7_3", len(chr3_cul_pos)))
#chr+=list(itertools.repeat("SM_V7_4", len(chr4_cul_pos)))
#chr+=list(itertools.repeat("SM_V7_5", len(chr5_cul_pos)))
#chr+=list(itertools.repeat("SM_V7_6", len(chr6_cul_pos)))
chr+=list(itertools.repeat("SM_V7_7", len(chr7_cul_pos)))


#ERR103048	bov
#ERR310937	curs
#ERR310940	marg

vcftools \
    --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
    --indv ERR103048 \
    --indv ERR310940 \
    --recode \
    --stdout \
    >bov_marg.vcf

 grep -v "#" bov_marg.vcf | cut -f3,10,11 | awk '{print $1"\t"substr($2,1,3)"\t"substr($3,1,3)}' | 
>bov_marg.gt

grep -P "1/1\t0/0" bov_marg.gt | cut -f1 >bov_marg.abba                                           
grep -P "0/0\t1/1" bov_marg.gt | cut -f1  >>bov_marg.abba  



vcftools \
    --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
    --snps bov_marg.abba \
    --keep niger.samples \
    --stdout \
    --freq \
    >niger.freq

vcftools \
    --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
    --snps bov_marg.abba \
    --keep tz.samples \
    --stdout \
    --freq \
    >tz.freq

vcftools \
    --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
    --snps bov_marg.abba \
    --indv ERR103048 \
    --stdout \
    --freq \
    >bov.freq


vcftools \
    --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
    --snps bov_marg.abba \
    --indv ERR310940 \
    --stdout \
    --freq \
    >marg.freq


paste niger.freq tz.freq bov.freq marg.freq | cut -f1,2,5,11,17,23 >freqs.out

#used excel to get a good table of abba baba sites: use them to extract from VCF



