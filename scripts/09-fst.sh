#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir fst

cd fst


#in python

import allel
import numpy as np  
import sys
import itertools

#read in the vcf data
callset=allel.read_vcf('../build_snp_panel/auto_maf.vcf', log=sys.stdout)
gt=allel.GenotypeArray(callset['calldata/GT'])

#get allele counts for each locus
ac=gt.count_alleles()


#designate the population (index) in the allele count array
egypt_pop=0 
bov_pop=1
mat_pop=[2,7,8]
guin_pop=3
inter_pop=4
curs_pop=5
marg_pop=6      
niger_pop=list(range(9,57))                                                            
tz_pop=list(range(58,103))

bov_group_pop=[bov_pop, curs_pop, marg_pop]
haem_group_pop=[egypt_pop + niger_pop + tz_pop]

#now count the alleles in the array for each pop
egypt_ac=gt.count_alleles(subpop=[egypt_pop])
bov_ac=gt.count_alleles(subpop=[bov_pop])
mat_ac=gt.count_alleles(subpop=mat_pop)
guin_ac=gt.count_alleles(subpop=[guin_pop])
inter_ac=gt.count_alleles(subpop=[inter_pop])
curs_ac=gt.count_alleles(subpop=[curs_pop])
marg_ac=gt.count_alleles(subpop=[marg_pop])
niger_ac=gt.count_alleles(subpop=niger_pop)
tz_ac=gt.count_alleles(subpop=tz_pop)

bov_group_ac=gt.count_alleles(subpop=bov_group_pop)
haem_group_ac=gt.count_alleles(subpop=haem_group_pop)

#calculate fst between populations
# via Willing et al. (2012) can't do FST with lt 4-6 samples
# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0042649

subpops=[niger_pop, tz_pop]
allel.average_weir_cockerham_fst(gt, subpops, blen=100)
#fst 0.2820221632974806
#se 0.010486528819494833
allel.average_patterson_fst(niger_ac, tz_ac, blen=100)
#fst 0.2857281133918937
#se 0.010557928500715526


#fst sliding window (niger vs. tz)----------------------------------------------
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

subpops=[niger_pop, tz_pop]
pos=callset['variants/POS']
chr1_fst_window, chr1_windows, chr1_counts=allel.windowed_weir_cockerham_fst(pos[chr1_start_idx:chr1_end_idx], gt[chr1_start_idx:chr1_end_idx], subpops, size=250000, step=50000, fill="Nan")   
chr2_fst_window, chr2_windows, chr2_counts=allel.windowed_weir_cockerham_fst(pos[chr2_start_idx:chr2_end_idx], gt[chr2_start_idx:chr2_end_idx], subpops, size=250000, step=50000, fill="Nan")   
chr3_fst_window, chr3_windows, chr3_counts=allel.windowed_weir_cockerham_fst(pos[chr3_start_idx:chr3_end_idx], gt[chr3_start_idx:chr3_end_idx], subpops, size=250000, step=50000, fill="Nan")   
chr4_fst_window, chr4_windows, chr4_counts=allel.windowed_weir_cockerham_fst(pos[chr4_start_idx:chr4_end_idx], gt[chr4_start_idx:chr4_end_idx], subpops, size=250000, step=50000, fill="Nan")   
chr5_fst_window, chr5_windows, chr5_counts=allel.windowed_weir_cockerham_fst(pos[chr5_start_idx:chr5_end_idx], gt[chr5_start_idx:chr5_end_idx], subpops, size=250000, step=50000, fill="Nan")   
chr6_fst_window, chr6_windows, chr6_counts=allel.windowed_weir_cockerham_fst(pos[chr6_start_idx:chr6_end_idx], gt[chr6_start_idx:chr6_end_idx], subpops, size=250000, step=50000, fill="Nan")   
chr7_fst_window, chr7_windows, chr7_counts=allel.windowed_weir_cockerham_fst(pos[chr7_start_idx:chr7_end_idx], gt[chr7_start_idx:chr7_end_idx], subpops, size=250000, step=50000, fill="Nan")   



#find average number of snps per window per chr
chr1_counts=list(chr1_counts)
chr2_counts=list(chr2_counts)
chr3_counts=list(chr3_counts)
chr4_counts=list(chr4_counts)
chr5_counts=list(chr5_counts)
chr6_counts=list(chr6_counts)
chr7_counts=list(chr7_counts)

chr1_avg_pos=list((chr1_windows[:,1] + chr1_windows[:,0])/2)
chr2_avg_pos=list((chr2_windows[:,1] + chr2_windows[:,0])/2)
chr3_avg_pos=list((chr3_windows[:,1] + chr3_windows[:,0])/2)
chr4_avg_pos=list((chr4_windows[:,1] + chr4_windows[:,0])/2)
chr5_avg_pos=list((chr5_windows[:,1] + chr5_windows[:,0])/2)
chr6_avg_pos=list((chr6_windows[:,1] + chr6_windows[:,0])/2)
chr7_avg_pos=list((chr7_windows[:,1] + chr7_windows[:,0])/2)

#numbers from the VCF ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf
chr1_len=88881357
chr2_len=48130368
chr3_len=50458499
chr4_len=47279781
chr5_len=25256119
chr6_len=24989083
chr7_len=19288021

chr1_cul_pos=[x+0 for x in chr1_avg_pos]
start=chr1_len
chr2_cul_pos=[x+start for x in chr2_avg_pos]
start+=chr2_len
chr3_cul_pos=[x+start for x in chr3_avg_pos]
start+=chr3_len
chr4_cul_pos=[x+start for x in chr4_avg_pos]
start+=chr4_len
chr5_cul_pos=[x+start for x in chr5_avg_pos]
start+=chr5_len
chr6_cul_pos=[x+start for x in chr6_avg_pos]
start+=chr6_len
chr7_cul_pos=[x+start for x in chr7_avg_pos]

cul_pos=chr1_cul_pos
cul_pos+=chr2_cul_pos
cul_pos+=chr3_cul_pos
cul_pos+=chr4_cul_pos
cul_pos+=chr5_cul_pos
cul_pos+=chr6_cul_pos
cul_pos+=chr7_cul_pos


chr=list(itertools.repeat("SM_V7_1", len(chr1_avg_pos)))
chr+=list(itertools.repeat("SM_V7_2", len(chr2_avg_pos)))
chr+=list(itertools.repeat("SM_V7_3", len(chr3_avg_pos)))
chr+=list(itertools.repeat("SM_V7_4", len(chr4_avg_pos)))
chr+=list(itertools.repeat("SM_V7_5", len(chr5_avg_pos)))
chr+=list(itertools.repeat("SM_V7_6", len(chr6_avg_pos)))
chr+=list(itertools.repeat("SM_V7_7", len(chr7_avg_pos)))


fst_per_window=list(chr1_fst_window) + list(chr2_fst_window) + list(chr3_fst_window) + list(chr4_fst_window) + list(chr5_fst_window) + list(chr6_fst_window) + list(chr7_fst_window) 
window_pos=chr1_avg_pos + chr2_avg_pos + chr3_avg_pos + chr4_avg_pos + chr5_avg_pos + chr6_avg_pos + chr7_avg_pos
counts_per_window=list(chr1_counts) + list(chr2_counts) + list(chr3_counts) + list(chr4_counts) + list(chr5_counts) + list(chr6_counts) + list(chr7_counts)

wc_fst=np.column_stack((chr, cul_pos, window_pos, fst_per_window, counts_per_window))

# remove if count is lt 10 #######

#save to file
np.savetxt("wc_fst_250kb_50kb_autosomal.csv", wc_fst, fmt='%s', delimiter=",")

#now plot this in R so that each chromosome has a different color

