#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# 17-diversity.sh - compare heterozygosity between tz and ne, and caclulate pi
#       in windows to build a ratio across chr4

# Uses a conda to manage the enviroment

#Set up the environment
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir diversity

cd diversity

cp ../niger.list .
cp ../tz.list .

#-------------------
#calculate heterozygosity
vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --maf 0.05 \
    --stdout \
    --het \
    --keep niger.list \
    --keep tz.list \
    >tz-ne_maf05.het

#plot in R

#-------------------
#calculate pi (diversity) on chr4

#get list of alleles with maf >0.05
vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --maf 0.05 \
    --keep niger.list \
    --keep tz.list \
    --stdout \
    --recode \
    >maf05.vcf

#create list of snps
grep -v "#" maf05.vcf | cut -f3 >maf05.list

#calculate pi in windows on chr 4 with snps from maf>0.05
#...in NE
vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --snps maf05.list \
    --keep niger.list \
    --chr SM_V7_4 \
    --window-pi 10000 \
    --window-pi-step 2500 \
    --stdout \
    >ne_SM_V7_4_maf05.window_pi

#...in TZ
vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --snps maf05.list \
    --keep tz.list \
    --chr SM_V7_4 \
    --window-pi 10000 \
    --window-pi-step 2500 \
    --stdout \
    >tz_SM_V7_4_maf05.window_pi

#combine the pi calculations to compare
paste ne_SM_V7_4_maf05.window_pi tz_SM_V7_4_maf05.window_pi >window.pi

#plot in R

