####################################################################################################

#to identify regions of introgression from an unknown source
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir sprime

cd sprime


grep "#" ../beagle_impute/cohort_snps_schMan_autosomal_beagle.vcf | tail -n1 | cut -f10- | sed 's/\t/\n/g' | grep "Sh.NE" >niger.list
grep "#" ../beagle_impute/cohort_snps_schMan_autosomal_beagle.vcf | tail -n1 | cut -f10- | sed 's/\t/\n/g' | grep "Sh.TZ" >tz.list


#only indentify introgressed regions (for now) in niger
vcftools \
    --vcf ../beagle_impute/cohort_snps_schMan_autosomal_beagle.vcf \
    --max-alleles 2 \
    --min-alleles 2 \
    --keep tz.list \
    --keep niger.list \
    --recode \
    --recode-INFO-all \
    --stdout >tz-niger_snps_schMan_autosomal_biallelic_beagle.vcf

grep -v "#" tz-niger_snps_schMan_autosomal_biallelic_beagle.vcf \
    | awk '{printf "%s\t%s\t%.6f\t%s\n", $1, $3, $2/287000, $2}' \
    >tz-niger_snps_schMan_autosomal_biallelic_beagle.map

#run sprime with niger vs. tz vcf and all tz as the outgroup
java -jar $WORK_DIR/scripts/sprime.11Apr18.bad.jar \
    gt=tz-niger_snps_schMan_autosomal_biallelic_beagle.vcf  \
    outgroup=tz.list \
    map=../beagle_impute/cohort_snps_schMan_autosomal_panel_single.map \
    out=niger_tz-out_sprime_2018-06-12 \
    mu=8.1E-9 \
    minscore=1.5e6


#exctract the introgressed snps into a bed file
#python build_sprime_intervals.py SM_V7_auto_sprime_2018-04-25.score >introgressed_segments.bed
#sed 1d SM_V7_auto_sprime_2018-04-25.score | awk '{print $3}' >introgressed_snps.list

#for speed - did this in excel
