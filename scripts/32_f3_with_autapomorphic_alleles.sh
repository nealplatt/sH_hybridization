source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir d_bovis_vs_curassoni

cd d_bovis_vs_curassoni

#D stats indicate introgression into nigerien haematobium from bovis and
# curassoni.  I am thinking this is driven by shared variation rather than
# introgression from both species.  Need to calculate D from autopomorphic loci.

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


#calcuate D with scikit allele. IN PYTHON !!!!!!!!!!!!!!!!!!!
import allel
import numpy as np  
import sys
import itertools

#read in the vcf data
callset=allel.read_vcf('bov_cur_autapomorphic_snps.vcf', log=sys.stdout)
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


allel.average_patterson_f3(niger_ac, tz_ac, bov_ac, 100)
#f3 -0.23441137635636838
#SE 0.012893599023468903
#Z  -18.18044565599514

allel.average_patterson_f3(niger_ac, tz_ac, curs_ac, 100)
#f3 -0.033345394285337684
#SE 0.024479469747819373
#Z  -1.362177964998939

allel.average_patterson_f3(tz_ac, niger_ac, bov_ac, 100)
#f3 3.164576701192813
#SE 0.44077365145752845
#Z  7.179595900817455

allel.average_patterson_f3(tz_ac, niger_ac, curs_ac, 100)
#f3 1.751804156402005
#SE 0.21958095642186146
#Z  7.977942099115457


#these results indicate that niger haem is an admixed tx and bov population.
# the test for admixture between tz and curs is negative (indicating admixture..
# but is not significant.
 
