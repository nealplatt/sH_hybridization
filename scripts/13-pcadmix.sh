#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# 13-pcadmix.sh - estimate snp ancestry with pcadmix 

# Uses a conda to manage the enviroment

#Set up the environment
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir pcadmix

cd pcadmix

#needs three seperate files (in bgl format)
#  - ancestral haplptype (bovis)
#  - ancestral haplptype (tz)
#  - admixed (NE)

#copy sample list files
cp ../tz.list .
cp ../niger.list .

#extract bovis in bgl
vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --indv ERR103048 \
    --recode \
    --stdout \
    | java -jar ../../scripts/vcf2beagle.jar NA BOV

#extract TZ in bgl
vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --keep tz.list \
    --recode \
    --stdout \
    | java -jar ../../scripts/vcf2beagle.jar NA TZ

#extract NE in bgl
vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --keep niger.list \
    --recode \
    --stdout \
    | java -jar ../../scripts/vcf2beagle.jar NA NE

#compress
gunzip BOV.bgl.gz TZ.bgl.gz NE.bgl.gz

#run pcadmix
$WORK_DIR/scripts/PCAdmix/PCAdmix3_linux \
    -anc TZ.bgl BOV.bgl \
    -adm NE.bgl \
    -o tz-bov-ne_maf00 \
    -lab TZ BOV NE \
    -w 30 \
    -thr 0.99 \
    -ld 0

#transform the output files into csv
sed 's/ /,/gi' tz-bov-ne_maf00.fbk.txt >tz-bov-ne_maf00.fbk.csv
sed 's/ /,/gi' tz-bov-ne_maf00.pca.txt >tz-bov-ne_maf00.pca.csv
sed 's/ /,/gi' tz-bov-ne_maf00.vit.txt >tz-bov-ne_maf00.vit.csv
sed 's/ /,/gi' tz-bov-ne_maf00.markers.txt >tz-bov-ne_maf00.markers.csv
sed 's/ /,/gi' tz-bov-ne_maf00.ia.txt >tz-bov-ne_maf00.ia.csv

############## with a few samples as a TZ control
# i wanted to test to see if we included some TZ as control...how much "bovis"
# ancestry would they have.  so i extracted 5 tz samples at random and then 
# included them in the admixed bgl file
#

shuf -n 5 tz.list >random.list
#Sh.TZ_PEM0125.1
#Sh.TZ_PEM0171.1
#Sh.TZ_PEM0103.1
#Sh.TZ_UNG0087.2
#Sh.TZ_UNG0076.1

#remove random samples
grep -v -f random.list tz.list >control.list

#add random samples
cat random.list niger.list >admixed.list

#extract CONTROL (tz minus random samples) and convert to bgl
vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --keep control.list \
    --recode \
    --stdout \
    | java -jar ../../scripts/vcf2beagle.jar NA control

#extract NE and random to make a "admixed" beagle file with controls
vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --keep admixed.list \
    --recode \
    --stdout \
    | java -jar ../../scripts/vcf2beagle.jar NA admixed

#compress
gunzip control.bgl.gz admixed.bgl.gz

#run pcadmix
$WORK_DIR/scripts/PCAdmix/PCAdmix3_linux \
    -anc control.bgl BOV.bgl \
    -adm admixed.bgl \
    -o control-bov-admixed_maf00 \
    -lab TZ BOV ADM \
    -w 30 \
    -thr 0.99 \
    -ld 0

#convert to csv
sed 's/ /,/gi' control-bov-admixed_maf00.fbk.txt >control-bov-admixed_maf00.fbk.csv
sed 's/ /,/gi' control-bov-admixed_maf00.pca.txt >control-bov-admixed_maf00.pca.csv
sed 's/ /,/gi' control-bov-admixed_maf00.vit.txt >control-bov-admixed_maf00.vit.csv
sed 's/ /,/gi' control-bov-admixed_maf00.markers.txt >control-bov-admixed_maf00.markers.csv
sed 's/ /,/gi' control-bov-admixed_maf00.ia.txt >control-bov-admixed_maf00.ia.csv

#now plot each snp in R
#x=chr and pos
#y=window_avg (avg freq in ne pop)
#color is based on window average
#id     chr     pos     window      window_score

