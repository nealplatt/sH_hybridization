#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# theta.sh - calculate theta with scikit-allele

# Uses a conda to manage environment

# this code is not intended to be run as is since it mixes bash and python
#   instead it reflects steps i used for the analysis

# Some steps run on local computer manually (format conversion)

#Set up the environment
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir theta

cd theta

#calcuate theta with scikit allele. IN PYTHON !!!!!!!!!!!!!!!!!!!
import allel
import numpy as np  
import sys
import itertools

#read in the vcf data
callset=allel.read_vcf('../beagle/auto_beagle_maf05.vcf', log=sys.stdout)
gt=allel.GenotypeArray(callset['calldata/GT'])

#get allele counts for each locus
ac=gt.count_alleles()

#designate the population (index) in the allele count array
egypt_pop=0 
bov_pop=1
curs_pop=2
niger_pop=list(range(3,48))                                                            
tz_pop=list(range(49,95))

#now count the alleles in the array for each pop
egypt_ac=gt.count_alleles(subpop=[egypt_pop])
bov_ac=gt.count_alleles(subpop=[bov_pop])
curs_ac=gt.count_alleles(subpop=[curs_pop])
niger_ac=gt.count_alleles(subpop=niger_pop)
tz_ac=gt.count_alleles(subpop=tz_pop)

#get ac for each chr [nums came from snps per chrom]
chr1_pos=list(range(1, 9501))
chr2_pos=list(range(9501, 9501+6135))
chr3_pos=list(range(9501+6135, 9501+6135+4592))
chr4_pos=list(range(9501+6135+4592, 9501+6135+4592+7677))
chr5_pos=list(range(9501+6135+4592+7677, 9501+6135+4592+7677+6473))
chr6_pos=list(range(9501+6135+4592+7677+6473, 9501+6135+4592+7677+6473+5744))
chr7_pos=list(range(9501+6135+4592+7677+6473+5744, 9501+6135+4592+7677+6473+5744+1559))

niger_theta=[]
niger_theta.append(allel.stats.watterson_theta(chr1_pos, niger_ac, start=chr1_pos[1], stop=chr1_pos[-1]))
niger_theta.append(allel.stats.watterson_theta(chr2_pos, niger_ac, start=chr2_pos[1], stop=chr2_pos[-1]))
niger_theta.append(allel.stats.watterson_theta(chr3_pos, niger_ac, start=chr3_pos[1], stop=chr3_pos[-1]))
niger_theta.append(allel.stats.watterson_theta(chr4_pos, niger_ac, start=chr4_pos[1], stop=chr4_pos[-1]))
niger_theta.append(allel.stats.watterson_theta(chr5_pos, niger_ac, start=chr5_pos[1], stop=chr5_pos[-1]))
niger_theta.append(allel.stats.watterson_theta(chr6_pos, niger_ac, start=chr6_pos[1], stop=chr6_pos[-1]))
niger_theta.append(allel.stats.watterson_theta(chr7_pos, niger_ac, start=chr7_pos[1], stop=chr7_pos[-1]))

niger_theta
#chr1 0.1846439591860867
#chr2 0.18647736891195757
#chr3 0.1863156319975177
#chr4 0.18377270540768453
#chr5 0.18603100421850133
#chr6 0.18605758253324417
#chr7 0.18199458728871643

#avg theta
theta=sum(niger_theta)/len(niger_theta)
#0.1850418342205298

#estimate Ne
mu=8.1e-9
ne=theta/(4*mu)
#ne = 5,711,167.7 or 5.71e6

tz_theta=[]
tz_theta.append(allel.stats.watterson_theta(chr1_pos, tz_ac, start=chr1_pos[1], stop=chr1_pos[-1]))
tz_theta.append(allel.stats.watterson_theta(chr2_pos, tz_ac, start=chr2_pos[1], stop=chr2_pos[-1]))
tz_theta.append(allel.stats.watterson_theta(chr3_pos, tz_ac, start=chr3_pos[1], stop=chr3_pos[-1]))
tz_theta.append(allel.stats.watterson_theta(chr4_pos, tz_ac, start=chr4_pos[1], stop=chr4_pos[-1]))
tz_theta.append(allel.stats.watterson_theta(chr5_pos, tz_ac, start=chr5_pos[1], stop=chr5_pos[-1]))
tz_theta.append(allel.stats.watterson_theta(chr6_pos, tz_ac, start=chr6_pos[1], stop=chr6_pos[-1]))
tz_theta.append(allel.stats.watterson_theta(chr7_pos, tz_ac, start=chr7_pos[1], stop=chr7_pos[-1]))

tz_theta
#chr1 0.10710213483233065
#chr2 0.09637082801869222
#chr3 0.10579642903276905
#chr4 0.10396907747439366
#chr5 0.09773849200607108
#chr6 0.10098346800502794
#chr7 0.11819905153174458


#avg theta
sum(tz_theta)/len(tz_theta)
#0.10430849727157561

theta=sum(tz_theta)/len(tz_theta)
#0.1850418342205298

#estimate Ne
mu=8.1e-9
ne=theta/(4*mu)
#Ne = 3,219,398.1 or 3.22e6
 
