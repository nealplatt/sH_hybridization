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

with open('f3_tz_ne_bov', 'w') as f:
    for item in f3_tz_ne_bov[3]:
        f.write("%s\n" % item)

with open('f3_ne_tz_bov', 'w') as f:
    for item in f3_ne_tz_bov[3]:
        f.write("%s\n" % item)

with open('f3_tz_ne_curs', 'w') as f:
    for item in f3_tz_ne_curs[3]:
        f.write("%s\n" % item)

with open('f3_ne_tz_curs', 'w') as f:
    for item in f3_ne_tz_curs[3]:
        f.write("%s\n" % item)




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



with open('af3_tz_ne_bov', 'w') as f:
    for item in af3_tz_ne_bov[3]:
        f.write("%s\n" % item)

with open('af3_ne_tz_bov', 'w') as f:
    for item in af3_ne_tz_bov[3]:
        f.write("%s\n" % item)

with open('af3_tz_ne_curs', 'w') as f:
    for item in af3_tz_ne_curs[3]:
        f.write("%s\n" % item)

with open('af3_ne_tz_curs', 'w') as f:
    for item in af3_ne_tz_curs[3]:
        f.write("%s\n" % item)


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

import allel
import numpy as np  
import sys
import itertools

#read in the vcf data

callset=allel.read_vcf('../build_snp_panel/auto_maf_ld.vcf', log=sys.stdout)
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
niger_pop=list(range(9,56))                                                            
tz_pop=list(range(56,102))

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


f3_tz_ne_bov=allel.average_patterson_f3(tz_ac, niger_ac, mat_ac, 50)

