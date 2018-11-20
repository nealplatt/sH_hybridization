####################################################################################################

#to identify regions of introgression from an unknown source
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir sprime

cd sprime

cp ../niger.list .
cp ../tz.list .

#only indentify introgressed regions (for now) in niger
vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --max-alleles 2 \
    --min-alleles 2 \
    --keep tz.list \
    --keep niger.list \
    --maf 0.031 \
    --recode \
    --recode-INFO-all \
    --stdout >tz-niger_auto_beagle.vcf

grep -v "#" tz-niger_auto_beagle.vcf \
    | awk '{printf "%s\t%s\t%.6f\t%s\n", $1, $3, $2/287000, $2}' \
    >tz-niger_auto_beagle.map

#run sprime with niger vs. tz vcf and all tz as the outgroup
java -jar $WORK_DIR/scripts/sprime.11Apr18.bad.jar \
    gt=tz-niger_auto_beagle.vcf  \
    outgroup=tz.list \
    map=tz-niger_auto_beagle.map \
    out=niger_tz-out_auto_sprime_2018-08-13 \
    mu=8.1E-9 \
    minscore=1e5


#exctract the introgressed snps into a bed file
python ../../scripts/build_sprime_intervals.py \
    niger_tz-out_auto_sprime_2018-08-13.score \
    >niger_tz-out_auto_sprime_2018-08-13.bed 

#create list of potentially introgressed snps
sed 1d niger_tz-out_auto_sprime_2018-08-13.score \
    | awk '{print $3}' \
    >niger_tz-out_auto_sprime_2018-08-13.list

#create list of SNPs id'd by sprime
cat niger_tz-out_auto_sprime_2018-08-13.score \
    | cut -f3 \
    | sed 1d \
    | sort \
    | uniq \
    >niger_tz-out_auto_sprime_2018-08-13.snps

#now do the reverse (how many novel segments in tanzanian samples)
java -jar $WORK_DIR/scripts/sprime.11Apr18.bad.jar \
    gt=tz-niger_auto_beagle.vcf  \
    outgroup=niger.list \
    map=tz-niger_auto_beagle.map \
    out=tz_niger-out_auto_sprime_2018-07-16 \
    mu=8.1E-9 \
    minscore=1e5


