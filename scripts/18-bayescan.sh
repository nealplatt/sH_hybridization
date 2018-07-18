source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir bayescan

cd bayescan

cp ../niger.list .
cp ../tz.list .

vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --keep niger.list \
    --keep tz.list \
    --stdout \
    --recode \
    >tz-ne.vcf
#After filtering, kept 95 out of 97 Individuals
#After filtering, kept 370,770 out of a possible 370,770 Sites


vcftools \
    --vcf tz-ne.vcf \
    --maf 0.05 \
    --stdout \
    --recode \
    --min-alleles 2 \
    --max-alleles 2 \
    >tz-ne_maf05_bi.vcf
#After filtering, kept 95 out of 95 Individuals
#After filtering, kept 36,308 out of a possible 370,770 Sites


#converted to bayescan format wtih pgd spider on local computer
#vcf ->pgd ->bayescan  (tz-ne_maf05_bi.bayescan)
