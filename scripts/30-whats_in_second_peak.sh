#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# 30-whats_in_the_second_peak.sh - From our bayescan results there is a second
#       minor peak showing selection...what gene is in the second peak 

# Uses a conda and singularity to manage the enviroment; needs relativley high
#       mem compute (used Titan - 125Gb)
source activate snp_calling
source /master/nplatt/schisto_hybridization/scripts/set_env.sh

cd $RESULTS_DIR

mkdir whats_in_second_peak

cd whats_in_second_peak

cp ../niger.list .
cp ../tz.list .

#second sig peak is from snps 28934, 28941, and 28942

grep -v "#" ../bayescan/tz-ne_maf05_bi.vcf | head -n 28934 | tail -n 1
#SM_V7_5 23065697        KL252047.1:1569 <...>
grep -v "#" ../bayescan/tz-ne_maf05_bi.vcf | head -n 28941 | tail -n 1
#SM_V7_5 25012178        KL251698.1:18653 <...>
grep -v "#" ../bayescan/tz-ne_maf05_bi.vcf | head -n 28942 | tail -n 1
#SM_V7_5 25019983        KL251698.1:30255 <...>

#find ancestry percentage in pcadmix
grep -e KL252047.1:1569 -e KL251698.1:18653 -e KL251698.1:30255 ../pcadmix/control-bov-admixed_maf00.vit.txt
#Window4717      KL252047.1:3349 KL252047.1:3316 KL252047.1:3312 KL252047.1:3275 KL252047.1:3225 KL252047.1:3203 KL252047.1:3173 KL252047.1:3157 KL252047.1:3121 KL252047.1:3038 KL252047.1:3034 KL252047.1:2964 KL252047.1:2835 KL252047.1:2760 KL252047.1:2745 KL252047.1:2722 KL252047.1:2670 KL252047.1:2602 KL252047.1:2405 KL252047.1:2374 KL252047.1:1760 KL252047.1:1730 KL252047.1:1697 KL252047.1:1689 KL252047.1:1628 KL252047.1:1569 KL252047.1:1537 KL252047.1:1532 KL252047.1:1511 KL252047.1:1469 
#Window4743      AMPZ01033274.1:106 KL253045.1:3671 KL253045.1:3115 KL251698.1:2109 KL251698.1:2136 KL251698.1:2199 KL253046.1:2126 KL253046.1:3152 KL251698.1:7884 KL251698.1:7953 KL251698.1:7960 KL251698.1:8011 KL251698.1:17475 KL251698.1:17514 KL251698.1:17521 KL251698.1:18565 KL251698.1:18653 KL251698.1:27498 KL251698.1:27542 KL251698.1:27560 KL251698.1:27579 KL251698.1:27599 KL251698.1:27636 KL251698.1:27729 KL251698.1:30250 KL251698.1:30255 KL251698.1:30257 KL251698.1:30261 KL251698.1:30270 KL251698.1:30282

#what is the ancestry percentage of these: calcuated manually
#Window4717	0.708333333
#Window4743	0.916666667

#what is xpehh look like in these regions?
#not significant (xpehh = 0.824474699716169) SM_V7_5:23065697
#not significant (xpehh = -1.297798788) SM_V7_5:25012178
#not significant (xpehh = -1.356006741) SM_V7_5:25019983
# in the general region xpehh doesn't exceed -1.3-1.5


#create bed file to see overlapping genes
echo -e "SM_V7_5\t23065695\t23065697\nSM_V7_5\t25012177\t25012179\nSM_V7_5\t25019982\t25019984" >chr5_targets.bed

bedtools intersect -wo -a chr5_targets.bed -b ../../data/Sm_v7.0.bed | cut -f1,2,7 | cut -f1 -d"." | sort | uniq
#SM_V7_5 23065695        Smp_169290
#SM_V7_5 25012177        Smp_081820
#SM_V7_5 25019982        Smp_081820

#These genes are: 
# http://www.genedb.org/gene/Smp_169290...a hypotheical protein
# http://www.genedb.org/gene/Smp_081820...a synaptic glycoprotein sc2

#not considering Smp_169290 any further, but may consider looking into Smp_081820 but xpehh wasn't sig.





