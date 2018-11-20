#...............................................................................
#...............................................................................
#PCA - needs LD pruned snps and to be in the plink ped format

#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

mkdir $RESULTS_DIR/06-pca
cd $RESULTS_DIR/06-pca

plink \
    --vcf $RESULTS_DIR/04-filter/cohort_snps_schMan_final_autosomal_LD.vcf \
    --pca \
    --allow-extra-chr \
    --out cohort_snps_schMan_final_autosomal_LD_pca


grep "#" $RESULTS_DIR/04-filter/cohort_snps_schMan_final_autosomal_LD.vcf  \
    | tail -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >cohort_snps_schMan_final_autosomal_LD.samples

vcftools \
    --vcf  $RESULTS_DIR/04-filter/cohort_snps_schMan_final_autosomal_LD.vcf \
    --remove-indv ERR103048 \
    --remove-indv ERR103051 \
    --remove-indv ERR119613 \
    --remove-indv ERR310937 \
    --remove-indv ERR310940 \
    --remove-indv ERR539855 \
    --remove-indv ERR539857 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_schMan_final_autosomal_LD_haemOnly.vcf

plink \
    --vcf cohort_snps_schMan_final_autosomal_LD_haemOnly.vcf \
    --pca \
    --allow-extra-chr \
    --out cohort_snps_schMan_final_autosomal_LD_haemOnly_pca


grep "#" cohort_snps_schMan_final_autosomal_LD_haemOnly.vcf  \
    | tail -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >cohort_snps_schMan_final_autosomal_LD_haemOnly.samples

#R code saved locally as mds.R


