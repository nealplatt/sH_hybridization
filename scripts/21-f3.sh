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



allel.average_patterson_f3(tz_ac, niger_ac, bov_ac, 100)
#f3 0.4517691002956402
#SE 0.03543316969051088
#Z  12.749892381674945

allel.average_patterson_f3(niger_ac, tz_ac, bov_ac, 100)
#f3 -0.054453718360619706
#SE 0.009507627944607596
#Z  -5.72737160918292

allel.average_patterson_f3(tz_ac, niger_ac, curs_ac, 100)
#f3 0.44908378098920204
#SE 0.03553567963537079
#Z  12.637545858056477

allel.average_patterson_f3(niger_ac, tz_ac, curs_ac, 100)
#f3 -0.05252711049980318
#SE 0.00995347006333908
#Z  -5.277266135884872


##################################
# in bash
plink \
    --vcf ../beagle/auto_beagle_maf05.vcf \
    --allow-extra-chr \
    --indep-pairwise 25 5 0.20 \
    --out auto_beagle_maf05

vcftools \
    --vcf ../beagle/auto_beagle_maf05.vcf \
    --exclude auto_beagle_maf05.prune.out \
    --recode \
    --recode-INFO-all \
    --stdout \
    >auto_beagle_maf05_LD.vcf

vcftools \
    --vcf auto_beagle_maf05_LD.vcf \
    --het \
    --stdout >auto_beagle_maf05_f.out

#F NE  = -0.001748478
#F TZ  =  0.3607874
#H0 NE =  0.2759724
#HO TZ =  0.1760968
