#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# 21-f3_and_derived_f3.sh - calculate genome wide fst in the scikit-allel 
#   then recalculate f3 with only derived alleles. 

# Uses a conda to manage the enviroment

# *** IMPORTANT
# So this script is not inteded to be "run" as much as it reflect steps that I
#   took in the analysis.  Its a combination of bash and python

#Set up the environment
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir f3

cd f3

#in python

import allel
import numpy as np  
import sys
import itertools

#read in the vcf data

callset=allel.read_vcf('auto_beagle_maf05_LD.vcf', log=sys.stdout)
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



f3_tz_ne_bov=allel.average_patterson_f3(tz_ac, niger_ac, bov_ac, 50)
#f3 0.4517691002956402
#SE 0.03543316969051088
#Z  12.749892381674945

f3_ne_tz_bov=allel.average_patterson_f3(niger_ac, tz_ac, bov_ac, 50)
#f3 -0.054453718360619706
#SE 0.009507627944607596
#Z  -5.72737160918292

f3_tz_ne_curs=allel.average_patterson_f3(tz_ac, niger_ac, curs_ac, 50)
#f3 0.44908378098920204
#SE 0.03553567963537079
#Z  12.637545858056477

f3_ne_tz_curs=allel.average_patterson_f3(niger_ac, tz_ac, curs_ac, 50)
#f3 -0.05252711049980318
#SE 0.00995347006333908
#Z  -5.277266135884872


################################################################################
mkdir f3_bovis_vs_curassoni

cd f3_bovis_vs_curassoni

#f3 stats indicates introgression from two populations in a third.  F3 is neg
# for NE, TZ + BOV and NE, TZ + CURS.  I am thinking this is driven by shared 
# variation rather than admixture from both species. 
#
# to test get autapomorphic/derived loci from each

#get list of autapomorphies
python ../../scripts/find_bovis_and_curassoni_autapomorphies.py \
    ../build_snp_panel/auto_maf.vcf \
    bov_cur_autapomorphic_snps.list

#get the vcf
cut -f1 bov_cur_autapomorphic_snps.list | sort | uniq >snps.list

vcftools \
    --vcf ../build_snp_panel/auto.vcf \
    --snps snps.list \
    --recode \
    --stdout \
    >bov_cur_autapomorphic_snps.vcf

###### in python
callset=allel.read_vcf('../f3_bovis_vs_curassoni/bov_cur_autapomorphic_snps.vcf', log=sys.stdout)
gt=allel.GenotypeArray(callset['calldata/GT'])

#get allele counts for each locus
ac=gt.count_alleles()

#designate the population (index) in the allele count array
bov_pop=1
curs_pop=5
marg_pop=6      
niger_pop=list(range(9,54))                                                            
tz_pop=list(range(55,101))

#now count the alleles in the array for each pop
bov_ac=gt.count_alleles(subpop=[bov_pop])
curs_ac=gt.count_alleles(subpop=[curs_pop])
marg_ac=gt.count_alleles(subpop=[marg_pop])
niger_ac=gt.count_alleles(subpop=niger_pop)
tz_ac=gt.count_alleles(subpop=tz_pop)


af3_ne_tz_bov=allel.average_patterson_f3(niger_ac, tz_ac, bov_ac, 100)
#f3 -0.23441137635636838
#SE 0.012893599023468903
#Z  -18.18044565599514

af3_ne_tz_curs=allel.average_patterson_f3(niger_ac, tz_ac, curs_ac, 100)
#f3 -0.033345394285337684
#SE 0.024479469747819373
#Z  -1.362177964998939

af3_tz_ne_bov=allel.average_patterson_f3(tz_ac, niger_ac, bov_ac, 100)
#f3 3.164576701192813
#SE 0.44077365145752845
#Z  7.179595900817455

af3_tz_ne_curs=allel.average_patterson_f3(tz_ac, niger_ac, curs_ac, 100)
#f3 1.751804156402005
#SE 0.21958095642186146
#Z  7.977942099115457


#these results indicate that niger haem is an admixed tx and bov population.
# the test for admixture between tz and curs is negative (indicating admixture..
# but is not significant.

