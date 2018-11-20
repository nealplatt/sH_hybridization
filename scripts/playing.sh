vcftools \
    --vcf cohort_snps_schMan_final_autosomal.vcf \
    --keep niger.list \
    --freq \
    --stdout \
    >niger.freq

vcftools \
    --vcf cohort_snps_schMan_final_autosomal.vcf \
    --keep zanzibar.list \
    --freq \
    --stdout \
    >zanzibar.freq

#margre
vcftools \
    --vcf cohort_snps_schMan_final_autosomal.vcf \
    --indv ERR310940 \
    --freq \
    --stdout \
    >margre.freq

#bovis
vcftools \
    --vcf cohort_snps_schMan_final_autosomal.vcf \
    --indv ERR103048 \
    --freq \
    --stdout \
    >bovis.freq

#curassoni
vcftools \
    --vcf cohort_snps_schMan_final_autosomal.vcf \
    --indv ERR310937 \
    --freq \
    --stdout \
    >curassoni.freq

 bcftools view \
    -m2 \
    -M2 \
    -v snps \
    cohort_snps_schMan_final_autosomal.vcf \
    >cohort_snps_schMan_final_autosomal_biallelic.vcf 


import allel
import numpy as np  
import sys

callset=allel.read_vcf('cohort_snps_schMan_final_autosomal_biallelic.vcf', fields=['numalt'], log=sys.stdout)
numalt=callset['variants/numalt']

np.max(numalt)                                                                      
#1
np.bincount(numalt)
array([     0, 346124])




callset=allel.read_vcf('cohort_snps_schMan_final_autosomal_biallelic.vcf', log=sys.stdout)
gt=allel.GenotypeArray(callset['calldata/GT'])

ac=gt.count_alleles()


marg_pop=6      
bov_pop=2 
niger_pop=list(range(10,58))                                                            
tz_pop=list(range(58,104))

marg_ac=gt.count_alleles(subpop=[marg_pop])
bov_ac=gt.count_alleles(subpop=[bov_pop])
niger_ac=gt.count_alleles(subpop=niger_pop)
tz_ac=gt.count_alleles(subpop=tz_pop)

chr1_end_idx=np.amax(np.where(callset['variants/CHROM']=="SM_V7_1"))
chr2_end_idx=np.amax(np.where(callset['variants/CHROM']=="SM_V7_2"))
chr3_end_idx=np.amax(np.where(callset['variants/CHROM']=="SM_V7_3"))
chr4_end_idx=np.amax(np.where(callset['variants/CHROM']=="SM_V7_4"))
chr5_end_idx=np.amax(np.where(callset['variants/CHROM']=="SM_V7_5"))
chr6_end_idx=np.amax(np.where(callset['variants/CHROM']=="SM_V7_6"))
chr7_end_idx=np.amax(np.where(callset['variants/CHROM']=="SM_V7_7"))

chr1_start_idx=np.amin(np.where(callset['variants/CHROM']=="SM_V7_1"))
chr2_start_idx=np.amin(np.where(callset['variants/CHROM']=="SM_V7_2"))
chr3_start_idx=np.amin(np.where(callset['variants/CHROM']=="SM_V7_3"))
chr4_start_idx=np.amin(np.where(callset['variants/CHROM']=="SM_V7_4"))
chr5_start_idx=np.amin(np.where(callset['variants/CHROM']=="SM_V7_5"))
chr6_start_idx=np.amin(np.where(callset['variants/CHROM']=="SM_V7_6"))
chr7_start_idx=np.amin(np.where(callset['variants/CHROM']=="SM_V7_7"))

chr1_idx=[chr1_start_idx, chr1_end_idx]
chr2_idx=[chr2_start_idx, chr2_end_idx]
chr3_idx=[chr3_start_idx, chr3_end_idx]
chr4_idx=[chr4_start_idx, chr4_end_idx]
chr5_idx=[chr5_start_idx, chr5_end_idx]
chr6_idx=[chr6_start_idx, chr6_end_idx]
chr7_idx=[chr7_start_idx, chr7_end_idx]

allel.average_weir_cockerham_fst(gt, subpops, blen, max_allele=None)

d=allel.patterson_d(tz_ac, niger_ac, bov_ac, marg_ac)

d_avg=allel.average_patterson_d(tz_ac, niger_ac, bov_ac, marg_ac, 10)

fst, se, block_fst, jackknife_fst=allel.average_weir_cockerham_fst(gt, subpop, 100)

allel.average_patterson_d(niger_ac, tz_ac, bov_ac, marg_ac, 100)

subpops=[niger_pop, tz_pop]

allel.moving_weir_cockerham_fst(gt, subpops, 100, step=50)
fst_window, windows, snp_counts=allel.windowed_weir_cockerham_fst(pos[chr1_start_idx:chr1_end_idx], gt[chr1_start_idx:chr1_end_idx], subpops, size=10000, step=5000, fill=nan)

fst_window, windows, snp_counts=allel.windowed_weir_cockerham_fst(pos[chr1_start_idx:chr1_end_idx], gt[chr1_start_idx:chr1_end_idx], subpops, size=500000, step=250000, fill="Nan")   
np.average(snp_counts)
#24.973574890376184
len(snp_counts)
#8666

avg_pos=(windows[:,1]+windows[:,0])/2

avg_pos=list(avg_pos)
fst_window=list(fst_window)

chr1_fst=np.column_stack((avg_pos, fst_window))

#save to file
np.savetxt("out.fst_5", chr1_fst, fmt='%s', delimiter=",")

#do for 

pop1=niger_pop
pop2=tz_pop

########################### FOLDED SFS #########################################




#outgroups to file
echo -e "ERR103051\nERR119613\nERR310940\nERR539855\nERR539857">outgroups.list
vcftools \
    --vcf cohort_snps_schMan_final_autosomal.vcf \
    --keep outgroups.list \
    --freq \
    --stdout \
    >outgroup.freq


cat outgroup.freq \
    | sed '1d' \
    | sed 's/:/\t/g' \
    | awk '{if ($8>$6) print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$5"\t"$6; \
        else \
        print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' \
    >tmp

awk '{if ($4>6 && $6==1) print $1"\t"$2}' tmp >anc.sites 
awk '{if ($4>6 && $6==1) print $1"\t"$2"\t"$5}' tmp >anc.alleles 

#need to find the ancestral allele.
add_anc_allele_info.py anc.alleles cohort_snps_schMan_final_autosomal.vcf aa.vcf

#get derived allele count
vcftools \
    --vcf aa.vcf \
    --derived \
    --keep niger.list \
    --counts2 \
    --stdout \
    >niger.counts2

vcftools \
    --vcf aa.vcf \
    --derived \
    --keep zanzibar.list \
    --counts2 \
    --stdout \
    >zanzibar.counts2

echo -e "chrom\tpos\tnigerCount\tzanzibarCount" >unfolded_sfs.out
paste niger.counts2 zanzibar.counts2 | sed '1d' | cut -f1,2,6,12 >>unfolded_sfs.out
cut -f3,4 unfolded_sfs.out | sed '1d' | sort | uniq -c | awk '{print $2","$3","$1}'>usfs_xy.csv


#ancestral allele appears first

#now using this information update the vcf with ancestral allele information in
# the info AA field


vcftools \
    --vcf cohort_snps_schMan_final_autosomal.vcf \
    --keep niger.list \
    --freq \
    --stdout \
    >niger.freq

vcftools \
    --vcf cohort_snps_schMan_final_autosomal.vcf \
    --keep zanzibar.list \
    --freq \
    --stdout \
    >zanzibar.freq

import allel
import numpy as np  
import sys

callset=allel.read_vcf('cohort_snps_schMan_final_autosomal_biallelic.vcf', log=sys.stdout)
gt=allel.GenotypeArray(callset['calldata/GT'])

d_callset=allel.read_vcf('aa.vcf', log=sys.stdout)
da=allel.GenotypeArray(callset['calldata/AA'])

ac=gt.count_alleles()


marg_pop=6      
bov_pop=2 
niger_pop=list(range(10,57))                                                            
tz_pop=list(range(58,104))

marg_ac=gt.count_alleles(subpop=[marg_pop])
bov_ac=gt.count_alleles(subpop=[bov_pop])
niger_ac=gt.count_alleles(subpop=niger_pop)
tz_ac=gt.count_alleles(subpop=tz_pop)


joint_fsfs=allel.joint_sfs_folded(tz_ac, niger_ac)

##############################################

paste zanzibar.freq niger.freq bovis.freq margre.freq     | cut -f1,2,5,11,17,23     | sed 's/.://g'     | awk '{if($6 == 1) {print $1"\t"$2"\t"1-$3"\t"1-$4"\t"1-$5"\t"1-$6} else print $0}'     | sed '1d'     | grep -vi nan     | awk '{if($6==0 && $5==1) print}' | awk '{if($3!=1 || $4!=1) print}' | awk '{if ($3!=0 || $4!=0
) print}' >bov_abba_babba                                                                  (snp_calling) [nplatt@compute-1-1633 play]$ paste zanzibar.freq niger.freq curassoni.freq margre.freq     | cut -f1,2,5,11,17,23     | sed 's/.://g'     | awk '{if($6 == 1) {print $1"\t"$2"\t"1-$3"\t"1-$4"\t"1-$5"\t"1-$6} else print $0}'     | sed '1d'     | grep -vi nan     | awk '{if($6==0 && $5==1) print}' | awk '{if($3!=1 || $4!=1) print}' | awk '{if ($3!=0 || $4!=0) print}' >cur_abba_babba       


plink \
    --vcf cohort_snps_schMan_final_autosomal.vcf \
    --out cohort_snps_schMan_final_autosomal \
    --transpose \
    --allow-extra-chr

plink \
    --vcf cohort_snps_schMan_final_autosomal.vcf \
    --out cohort_snps_schMan_final_autosomal \
    --allow-extra-chr \
    --recode

python2 popstats.py \
    --file cohort_snps_schMan_final_autosomal \
    --pops Sh.NE,Sh.TZ,ERR103048,ERR310940 \
    --informative \
    --not23

#outgroup
vcftools \
    --vcf cohort_snps_schMan_final_autosomal_LD.vcf \
    --indv ERR310940 \
    --freq \
    --stdout \
    >margre.freq
ERR310940


#create anc.fasta
angsd -doFasta 1 \
    -i ../processed_reads/sh_mapped_reads/ERR310940_processed.bam \
    -doCounts 1 \
    -out ERR310940

angsd -doFasta 1 \
    -i ../processed_reads/sh_mapped_reads/Sh.NE_Dai-010.1_processed.bam \
    -doCounts 1 \
    -out Sh.NE_Dai-010.1

gunzip ERR310940.fa.gz
samtools ERR310940.fa

#create file list of 
#NE, TZ, Bov and Curs
ls ../processed_reads/sh_mapped_reads/Sh.NE*bam >NE_bams.list
ls ../processed_reads/sh_mapped_reads/Sh.TZ*bam >TZ_bams.list
ls ../processed_reads/sh_mapped_reads/ERR103048*bam >BOV_bams.list
ls ../processed_reads/sh_mapped_reads/ERR310937*bam >CUR_bams.list

angsd -doFasta 1 \
    -i ../processed_reads/sh_mapped_reads/Sh.NE_Dai-010.1_processed.bam \
    -doCounts 1 \
    -out Sh.NE_Dai-010.1

gunzip Sh.NE_Dai-010.1.fa.gz
samtools Sh.NE_Dai-010.1.fa


angsd -doAbbababa2



python ABBABABAwindows.py \
    -g /zoo/disk1/shm45/vcf/set62/set62.chr21.DP5GQ30.AN100MAC1.diplo.gz \
    -f phased \
    -o output.csv \
    -w 100000 \
    -m 100 \
    -s 100000 \
    -p P1 A1,A2,A3,A4 \
    -p P2 B1,B2,B3,B4 \
    -p P3 C1,C2,C3,C4 \
    -p O D1,D2,D3,D4 \
    -T 10 \
    --minData 0.5

python parseVCF.py \
    -i cohort_snps_schMan_final.vcf \
    | gzip \
        > output.geno.gz

python ABBABABAwindows.py \
    -g /zoo/disk1/shm45/vcf/set62/set62.chr21.DP5GQ30.AN100MAC1.diplo.gz \
    -f phased \
    -o output.csv \
    -w 100000 \
    -m 100 \
    -s 100000 \
    -p P1 A1,A2,A3,A4 \
    -p P2 B1,B2,B3,B4 \
    -p P3 C1,C2,C3,C4 \
    -p O D1,D2,D3,D4 \
    -T 10 \
    --minData 0.5


import allel
import sys

callset=allel.read_vcf('cohort_snps_schMan_final_autosomal.vcf', log=sys.stdout)



g = allel.HaplotypeArray([[0, 0, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 1],
    [0, 1, 1, 1],
    [1, 1, 1, 1],
    [0, 0, 1, 2],
    [0, 1, 1, 2],
    [0, 1, -1, -1]])




bcftools filter \
    -m2 --threads 10 \
    cohort_snps_schMan_final.vcf \
    >polymorphic.vcf
-m2 -M2 -v snps
################## ARE THESE SITES POLYMORPHIC


#get autosomal data (from filtered data)
vcftools \
    --vcf cohort_snps_schMan_final.vcf \
    --chr SM_V7_1 \
    --chr SM_V7_2 \
    --chr SM_V7_3 \
    --chr SM_V7_4 \
    --chr SM_V7_5 \
    --chr SM_V7_6 \
    --chr SM_V7_7 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_schMan_final_autosomal.vcf
    #After filtering, kept 105 out of 105 Individuals
    #After filtering, kept 348131 out of a possible 458403 Sites


#get autosomal LD data
plink \
    --vcf cohort_snps_schMan_final_autosomal.vcf \
    --allow-extra-chr \
    --indep-pairwise 25 5 0.2\
    --out cohort_snps_schMan_final_autosomal_LD-25-5-2

vcftools \
    --vcf cohort_snps_schMan_final_autosomal.vcf \
    --snps cohort_snps_schMan_final_autosomal_LD-25-5-2.prune.in \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_schMan_final_autosomal_LD.vcf
    #After filtering, kept 105 out of 105 Individuals
    #After filtering, kept 102863 out of a possible 348131 Sites



#get mitochondrial data (from initial snp data --no site/indiv filtering)
vcftools \
    --vcf cohort_snps_annotated.vcf \
    --chr AMPZ01026399.1 \
    --recode \
    --recod-INFO-all \
    --stdout \
    >cohort_snps_mito.vcf
#After filtering, kept 117 out of 117 Individuals
#After filtering, kept 450 out of a possible 611210 Sites


vcftools \
    --vcf cohort_snps_schMan_final_autosomal_LD.vcf \
    --weir-fst-pop zanzibar.list \
    --weir-fst-pop niger.list







#FST
