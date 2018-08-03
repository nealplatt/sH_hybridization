#...............................................................................
#...............................................................................
#PCA - needs LD pruned snps and to be in the plink ped format

#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

mkdir $RESULTS_DIR/pca
cd $RESULTS_DIR/pca

#the goal is to end up with 3 runs of PCA
#1) all samples
#2) curs, bovis, haem
#3) haem only

#cohort_snps_schMan_maf-05-perc_autosomal_panel_LD-25-5-2.vcf

################################################################################
#1) ALL SAMPLES
plink \
    --vcf ../build_snp_panel/auto_maf_ld.vcf \
    --pca \
    --allow-extra-chr \
    --out auto_maf_ld_pca
    #After filtering, kept 104 out of 104 Individuals
    #After filtering, kept 5905 out of a possible 5905 Sites

#redo the pca but remove all non-haem/bovis/curassoni samples
grep "#" ../build_snp_panel/auto_maf_ld.vcf  \
    | tail -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >auto_maf_ld.samples

################################################################################
#2) HAEM_GROUP (includes bovis)

vcftools \
    --vcf  ../build_snp_panel/auto_maf_ld.vcf \
    --remove-indv ERR103051 \
    --remove-indv ERR119612 \
    --remove-indv ERR119613 \
    --remove-indv ERR310940 \
    --remove-indv ERR539855 \
    --remove-indv ERR539857 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >haem_auto_maf_ld.vcf
    #After filtering, kept 96 out of 104 Individuals
    #After filtering, kept 5905 out of a possible 5905 Sites

plink \
    --vcf haem_auto_maf_ld.vcf \
    --pca \
    --allow-extra-chr \
    --out haem_auto_maf_ld_pca


grep "#" haem_auto_maf_ld.vcf \
    | tail -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >haem_auto_maf_ld_pca.samples

################################################################################
#3) HAEM_ONLY (includes haem from SRA)
vcftools \
    --vcf haem_auto_maf_ld.vcf \
    --remove-indv ERR103048 \
    --remove-indv ERR310937 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >haemOnly_auto_maf_ld.vcf
    #After filtering, kept 94 out of 96 Individuals
    #After filtering, kept 5882 out of a possible 5882 Sites

plink \
    --vcf haemOnly_auto_maf_ld.vcf \
    --pca \
    --allow-extra-chr \
    --out haemOnly_auto_maf_ld_pca


grep "#" haemOnly_auto_maf_ld.vcf  \
    | tail -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >haemOnly_auto_maf_ld_pca.samples

# PLOT IN R (pca.R)
#R $SCRIPTS_DIR/R/pca.R
