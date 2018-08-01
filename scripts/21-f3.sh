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


allel.average_patterson_f3(tz_ac, niger_ac, bov_ac, 100)
#f3 1.3067408628912387
#SE 0.07403752556207435
#Z  17.6497101026917

allel.average_patterson_f3(niger_ac, tz_ac, bov_ac, 100)
#f3 -0.1629027774840651
#SE 0.007224948856017552
#Z  -22.54725683606546
