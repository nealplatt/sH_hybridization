#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir beagle

cd beagle

#calculate a recomb map.  1cm = 287Kb based on data from Criscione et al. (2009) avg
grep -v "#" ../build_snp_panel/auto.vcf  \
    | awk '{printf "%s\t%s\t%.6f\t%s\n", $1, $3, $2/287000, $2}' \
    >auto.map

#create a list of sample to imput genotypes - in this analysis keeping all
#  haem samples and then the "good" bovis sample
grep "#" ../build_snp_panel/auto.vcf \
    | tail -n 1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >samples.list

#nano samples.list (and remove unnecessary...outgroup...samples)

#impute on each seperatley
for i in $(seq 1 7); do

    CHR=SM_V7_"$i"

    #extract autosome specific vcf for the samples of interest
    vcftools \
        --vcf ../build_snp_panel/auto.vcf \
        --chr $CHR \
        --keep samples.list \
        --recode \
        --stdout >$CHR".vcf"

    sed -i 's/,assembly=schMan_v7.fa//gi' $CHR.vcf

    grep $CHR auto.map >$CHR".map"

    beagle \
        gt=$CHR".vcf" \
        out=$CHR"_beagle" \
        map=$CHR".map" \
        nthreads=4 \
        window=300 \
        overlap=30 \
        niterations=250 >>$CHR.beaglelog 2>&1 &

done

wait 

gunzip SM*vcf.gz

#create one file of phased snps (for all autosomes)
vcfcombine \
    SM_V7_1_beagle.vcf \
    SM_V7_2_beagle.vcf \
    SM_V7_3_beagle.vcf \
    SM_V7_4_beagle.vcf \
    SM_V7_5_beagle.vcf \
    SM_V7_6_beagle.vcf \
    SM_V7_7_beagle.vcf \
    >auto_beagle.vcf
    #370,770

vcftools \
    --vcf auto_beagle.vcf \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >auto_beagle_maf05.vcf
    #41,680
