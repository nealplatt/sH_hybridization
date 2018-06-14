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
vcftools \
    --vcf ../beagle_impute/cohort_snps_schMan_autosomal_beagle.vcf \
    --indv ERR103048 \
    --recode \
    --stdout \
    | java -jar ../../scripts/vcf2beagle.jar NA bovis_autosomal

#extract TZ
cat ../beagle_impute/cohort_snps_schMan_autosomal_beagle.vcf \
    | grep "#" \
    | tail -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    | grep Sh.TZ \
    >tz.samples

#extract NE
vcftools \
    --vcf ../beagle_impute/cohort_snps_schMan_autosomal_beagle.vcf \
    --keep tz.samples \
    --recode \
    --stdout \
    | java -jar ../../scripts/vcf2beagle.jar NA TZ_autosomal

#create a file of ONLY admixed individualsadmixed set
cat ../beagle_impute/cohort_snps_schMan_autosomal_beagle.vcf \
    | grep "#" \
    | tail -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    | grep Sh.NE \
    >ne.samples

vcftools \
    --vcf ../beagle_impute/cohort_snps_schMan_autosomal_beagle.vcf \
    --keep ne.samples \
    --recode \
    --stdout \
    | java -jar ../../scripts/vcf2beagle.jar NA NE_autosomal

gunzip *.gz

# NEED TO FIGURE OUT MAPS
##get the genetic map and physical map
plink \
    --vcf ../beagle_impute/cohort_snps_schMan_autosomal_beagle.vcf \
    --out cohort_snps_schMan_autosomal_beagle \
    --recode12 \
    --allow-extra-chr

cut -f4 cohort_snps_schMan_autosomal_beagle.map \
    | awk '{printf "%s\t-\t%.10f\n", $1, 287000/1000000}' \
    | head  \
    >cohort_snps_schMan_autosomal_beagle.geneticmap

$WORK_DIR/scripts/PCAdmix/PCAdmix3_linux \
    -anc TZ_autosomal.bgl bovis_autosomal.bgl \
    -adm NE_autosomal.bgl \
    -o num_2 \
    -bed num_2.bed \
    -lab TZ BOV NE_ADMIXED

$WORK_DIR/scripts/PCAdmix/PCAdmix3_linux \
    -anc TZ_autosomal.bgl bovis_autosomal.bgl \
    -adm NE_autosomal.bgl \
    -o num_2 \
    -bed num_2.bed \
    -lab TZ BOV NE_ADMIXED \
    -wSNP



shuf tz.samples | shuf | shuf | shuf | shuf | shuf | shuf | shuf  >tz_test.sampl 
head -n 20 tz_shuffled.samples >tz_test.samples
tail -n 27 tz_shuffled.samples >tz_control.samples
cat tz_test.samples ne.samples >adm.samples


vcftools \
    --vcf ../beagle_impute/cohort_snps_schMan_autosomal_beagle.vcf \
    --keep adm.samples \
    --chr SM_V7_4 \
    --recode \
    --stdout \
    | java -jar ../../scripts/vcf2beagle.jar NA ADM_autosomal

vcftools \
    --vcf ../beagle_impute/cohort_snps_schMan_autosomal_beagle.vcf \
    --keep tz_control.samples \
    --chr SM_V7_4 \
    --recode \
    --stdout \
    | java -jar ../../scripts/vcf2beagle.jar NA TZ_CONTROL_autosomal

gunzip *.gz

$WORK_DIR/scripts/PCAdmix/PCAdmix3_linux \
    -anc bovis_autosomal.bgl TZ_CONTROL_autosomal.bgl \
    -adm ADM_autosomal.bgl \
    -o test_adm \
    -map cohort_snps_schMan_autosomal_beagle.map \
    -bed \
    -lab BOVIS TZ NE \
    -wMB 0.25


$WORK_DIR/scripts/PCAdmix/PCAdmix3_linux


