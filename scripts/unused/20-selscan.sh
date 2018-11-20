source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir selscan

cd selscan

#create a list of NE and TZ samples
cp ../niger.list .
cp ../tz.list .
cat niger.list tz.list >samples.list


#convert to phased, binary (01), chromosome vcfs
for i in $(seq 1 7); do

    CHR=SM_V7_"$i"

    #extract autosome specific vcf for the samples of interest
    vcftools \
        --vcf ../beagle/auto_beagle_maf05.vcf \
        --chr $CHR \
        --keep samples.list \
        --recode \
        --stdout >$CHR".vcf"

    plink \
        --vcf $CHR".vcf" \
        --out $CHR \
        --recode12 \
        --allow-extra-chr

done


#calculate ihs
../../scripts/selscan/selscan \
    -ihs


#calculate xpEHH 
#need to create pop specific vcf files
for i in $(seq 1 7); do
    for POP in niger tz; do

        CHR=SM_V7_"$i"

        #extract pop specific files per chr
        vcftools \
            --vcf $CHR".vcf" \
            --keep $POP.list \
            --recode \
            --stdout >$CHR"_"$POP".vcf"
    done
done


../../scripts/selscan/selscan \
    --threads 10 \
    --ihs \
    --vcf SM_V7_1.vcf \
    --map SM_V7_1.map \
    --out SM_V7_1_ihs_def.test



    --cutoff
    --maf
    --gap-scale 
    --max-gap
--max-extend
--trunc-ok
--alt
--skip-low-freq
