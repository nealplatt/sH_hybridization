#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir beagle_impute

cd beagle_impute

#imputing snps at all autosomal loci (regardless of LD)
grep -v "#" ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
    | cut -f1,2 \
    | sort \
    | uniq -c \
    | awk '{ if($1>1) print $2"\t"$3}' \
    >remove.pos

#for some reason there are still a handful of duplicate snps positions - remove:
vcftools \
    --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
    --exclude-positions remove.pos \
    --recode \
    --recode-INFO-all \
    --stdout \
    | sed 's/,assembly=schMan_v7.fa//gi' \
    >cohort_snps_schMan_autosomal_panel_single.vcf
    #After filtering, kept 338881 out of a possible 338881 Sites

#calculate a recomb map.  1cm = 287Kb based on data from Criscione et al. (2009) avg
grep -v "#" cohort_snps_schMan_autosomal_panel_single.vcf \
    | awk '{printf "%s\t%s\t%.6f\t%s\n", $1, $3, $2/287000, $2}' \
    >cohort_snps_schMan_autosomal_panel_single.map

#create a list of sample to imput genotypes - in this analysis keeping all
#  haem samples and then the "good" bovis sample
grep "#" cohort_snps_schMan_autosomal_panel_single.vcf \
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
        --vcf cohort_snps_schMan_autosomal_panel_single.vcf \
        --chr $CHR \
        --keep samples.list \
        --recode \
        --stdout >$CHR.vcf

    sed -i 's/,assembly=schMan_v7.fa//gi' $CHR.vcf

    grep $CHR cohort_snps_schMan_autosomal_panel_single.map >$CHR.map

    beagle \
        gt=$CHR.vcf \
        out=$CHR"_beagle" \
        map=$CHR.map \
        nthreads=10 \
        lowmem=true \
        window=302 \
        overlap=30 \
        niterations=250

done

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
    >cohort_snps_schMan_autosomal_beagle.vcf



