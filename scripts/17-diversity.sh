#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir diversity

cd diversity

cp ../niger.list .
cp ../tz.list .

############################## R O H ###########################################
#identify regions of homozygosity in NE and TZ (and find regions that differ)
for i in $(seq 1 7); do
    CHR=SM_V7_"$i"
    for POP in niger tz; do
    
        vcftools \
            --vcf ../beagle/auto_beagle.vcf \
            --chr $CHR \
            --LROH \
            --keep $POP.list \
            --stdout \
            >$CHR"_"$POP.lroh

        sed '1d' $CHR"_"$POP.lroh >$CHR"_"$POP.bed
    done
done

#gen coords for plotting in R
cat SM_V7_?_tz.bed | awk '{print $5-$4","$6}' >tz.cords
cat SM_V7_?_ne.bed | awk '{print $5-$4","$6}' >ne.cords
cat ne.cords tz.cords >all.coords

############################## HET & F #########################################
vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --maf 0.05 \
    --stdout \
    --het \
    --keep niger.list \
    --keep tz.list \
    >tz-ne_maf05.het

#plot in R

