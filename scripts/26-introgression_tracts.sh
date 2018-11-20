#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

cp ../build_snp_panel/auto_maf.vcf .


#get tz vcf
vcftools \
    --vcf auto_maf.vcf \
    --keep ../tz.list \
    --freq \
    --stdout \
    >tz.freq

#get bovis vcf
vcftools \
    --vcf auto_maf.vcf \
    --indv ERR103048 \
    --freq \
    --stdout \
    >bov.freq

#get outgroup vcf
vcftools \
    --vcf auto_maf.vcf \
    --indv ERR103051 \
    --indv ERR119612 \
    --indv ERR119613 \
    --indv ERR310937 \
    --indv ERR310940 \
    --indv ERR539855 \
    --indv ERR539857 \
    --freq \
    --stdout \
    >outgroup.freq

paste bov.freq tz.freq outgroup.freq >freq
                                                                                                                                                                            


#find private bovis alleles (must remove NE samples)

vcftools \
    --vcf auto_maf.vcf \
    --remove ../niger.list \
    --remove-indv ERR037800 \
    --stdout \
    --recode \
    --min-alleles 2 \
    --max-alleles 2 \
    >tz_and_outs_het.vcf

#find private bovis alleles
 grep "#" | tail -n 1




#find private tz alleles (anything that is in tz but not in other no NE SH)
