#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir pcadmix

cd pcadmix

#needs admixed file
#and one file per ancestral halpotype

#create files in beagle format
#exctract bovis
cp ../tz.list .
cp ../niger.list .


vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --indv ERR103048 \
    --recode \
    --stdout \
    | java -jar ../../scripts/vcf2beagle.jar NA BOV

#extract TZ
vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --keep tz.list \
    --recode \
    --stdout \
    | java -jar ../../scripts/vcf2beagle.jar NA TZ

    #extract NE
vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --keep niger.list \
    --recode \
    --stdout \
    | java -jar ../../scripts/vcf2beagle.jar NA NE

gunzip BOV.bgl.gz TZ.bgl.gz NE.bgl.gz

$WORK_DIR/scripts/PCAdmix/PCAdmix3_linux \
    -anc TZ.bgl BOV.bgl \
    -adm NE.bgl \
    -o tz-bov-ne_maf00 \
    -lab TZ BOV NE \
    -w 30 \
    -thr 0.99 \
    -ld 0

sed 's/ /,/gi' tz-bov-ne_maf00.fbk.txt >tz-bov-ne_maf00.fbk.csv
sed 's/ /,/gi' tz-bov-ne_maf00.pca.txt >tz-bov-ne_maf00.pca.csv
sed 's/ /,/gi' tz-bov-ne_maf00.vit.txt >tz-bov-ne_maf00.vit.csv
sed 's/ /,/gi' tz-bov-ne_maf00.markers.txt >tz-bov-ne_maf00.markers.csv
sed 's/ /,/gi' tz-bov-ne_maf00.ia.txt >tz-bov-ne_maf00.ia.csv


############## with a few samples as a TZ control

shuf -n 5 tz.list
#Sh.TZ_PEM0125.1
#Sh.TZ_PEM0171.1
#Sh.TZ_PEM0103.1
#Sh.TZ_UNG0087.2
#Sh.TZ_UNG0076.1

#remove random samples
nano tz.list >control.list

#add random samples
nano niger.list >admixed.list

#extract CONTROL
vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --keep control.list \
    --recode \
    --stdout \
    | java -jar ../../scripts/vcf2beagle.jar NA control

    #extract NE
vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --keep admixed.list \
    --recode \
    --stdout \
    | java -jar ../../scripts/vcf2beagle.jar NA admixed

gunzip control.bgl.gz admixed.bgl.gz

$WORK_DIR/scripts/PCAdmix/PCAdmix3_linux \
    -anc control.bgl BOV.bgl \
    -adm admixed.bgl \
    -o control-bov-admixed_maf00 \
    -lab TZ BOV ADM \
    -w 30 \
    -thr 0.99 \
    -ld 0

sed 's/ /,/gi' control-bov-admixed_maf00.fbk.txt >control-bov-admixed_maf00.fbk.csv
sed 's/ /,/gi' control-bov-admixed_maf00.pca.txt >control-bov-admixed_maf00.pca.csv
sed 's/ /,/gi' control-bov-admixed_maf00.vit.txt >control-bov-admixed_maf00.vit.csv
sed 's/ /,/gi' control-bov-admixed_maf00.markers.txt >control-bov-admixed_maf00.markers.csv
sed 's/ /,/gi' control-bov-admixed_maf00.ia.txt >control-bov-admixed_maf00.ia.csv

# NEED TO FIGURE OUT MAPS
##get the genetic map and physical map
#plink \
#    --vcf ../beagle/auto_beagle.vcf \
#    --out auto_beagle \
#    --recode12 \
#    --allow-extra-chr

#cut -f4 auto_beagle.map \
#    | awk '{printf "%s\t-\t%.10f\n", $1, 287000/1000000}' \
#    | head  \
#    >auto_beagle.geneticmap



