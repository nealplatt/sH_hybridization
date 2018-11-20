#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# 22-prep_for_rehh-xpehh.sh - prep the vcf files into a format that is
#   compatible wtih rehh in R  

# Uses a conda to manage the enviroment

#Set up the environment
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir rehh

cd rehh

#create a list of NE and TZ samples
cp ../niger.list .
cp ../tz.list .
cat niger.list tz.list >samples.list


#now format the data for rehh
#change chr names from SM_V7_X to X
sed -s 's/SM_V7_//gi' ../beagle/auto_beagle.vcf >auto_beagle_chr.vcf

#creat inp and map files for each chromosome
for CHR in $(seq 1 7); do
    for POP in  tz niger; do
        
        vcftools \
            --vcf auto_beagle_chr.vcf \
            --keep $POP.list \
            --chr $CHR \
            --recode \
            --stdout >$POP"_SM_V7_"$CHR.vcf

        plink \
            --vcf $POP"_SM_V7_"$CHR.vcf \
            --out $POP \
            --recode fastphase 

        grep -v "#" $POP.chr-$CHR.recode.phase.inp \
            | sed -e 's/\(.\)/\1 /g' \
            | sed 1,3d \
            | awk '{print NR" "$0}' \
            >$POP"_"$CHR.inp
    
        grep -v "#" $POP"_SM_V7_"$CHR.vcf \
            | awk '{print $3"\t"$1"\t"$2"\t"$4"\t"$5}' \
            >$POP"_"$CHR.map
    done
done 

#combine into files suitable for whole genome analysis
for CHR in $(seq 1 7); do
    for POP in  tz niger; do
        
       grep -v "#" $POP.chr-$CHR.recode.phase.inp \
            | sed -e 's/\(.\)/\1 /g' \
            | sed 1,3d \
            >$POP"_"$CHR.tmp
    
    done
done 
 
#create genome wide inp files for each pop
paste niger_?.tmp | awk '{print NR" "$0}' >niger.inp
paste tz_?.tmp | awk '{print NR" "$0}' >tz.inp

#create genome wide map files for each pop
cat niger_?.map >niger.map
cat tz_?.map >tz.map

#now analyze in rehh in R

