# split cohort_r3_filteredVariants.g.vcf into indels and snps

source /master/nplatt/sH_hybridization/scripts/set-env.sh

rm m1*.vcf* m2*.vcf*
rm tmp_cohort.vcf*
rm *.list

mkdir contig_vcfs
mv AMP*.vcf* contig_vcfs 
mv KL*.vcf* contig_vcfs

mkdir individual_vcfs
mv Sh.*final.g.vcf* individual_vcfs/

cd /master/nplatt/sH_hybridization/results/filter_cohort_vcf

IN_RAW_VCF=/master/nplatt/sH_hybridization/results/haplotype_caller/cohort_raw_postBQSR.vcf

SNP_VCF=cohort_raw_SNPs_r3.vcf
INDEL_VCF=cohort_raw_INDELs_r3.vcf
REFERENCE=/master/nplatt/sH_hybridization/data/genome/schHae_v1.fa

FILTERED_INDELS_VCF=cohort_filtered_INDELs_r3.vcf
FILTERED_SNPS_VCF=cohort_filtered_SNPs_r3.vcf

MERGED_VARIANTS_VCF=cohort_filtered_variants_r3.vcf

#split vcfs into indels and snps
$SINGULARITY gatk SplitVcfs \
      --INPUT=$IN_RAW_VCF \
      --SNP_OUTPUT=$SNP_VCF \
      --INDEL_OUTPUT=$INDEL_VCF \
      --STRICT=false

#i am thinking it may be worth removing
#1) multiallelic sites
#2) sites that are NOT snps or indels
#3) should we remove low coverage samples
#4) frequent SNPS in dense regions (possibly pseudo gene)
#5) SNPs within X bp of an indel --mask --maskExtensions --maskName
#6) max no call samples per site

#filter SNPs
$SINGULARITY gatk VariantFiltration \
   -R $REFERENCE  \
   -O $FILTERED_SNPS_VCF \
   -V $SNP_VCF \
   --filter-name "snp_QD_lt_5" \
   --filter-expression "QD < 5.0" \
   --filter-name "snp_FS_gt_55" \
   --filter-expression "FS > 55.0" \
   --filter-name "snp_MQ_lt_40" \
   --filter-expression "MQ < 40.0" \
   --filter-name "snp_MQRankSum_lt_-12.5" \
   --filter-expression "MQRankSum < -12.5" \
   --filter-name "snp_ReadPosRankSum_lt_-8" \
   --filter-expression "ReadPosRankSum < -8.0" \
   --filter-name "snp_SQR_gt_3" \
   --filter-expression "SOR > 3.0" &

#filter INDELs
$SINGULARITY gatk VariantFiltration \
   -R $REFERENCE  \
   -V $INDEL_VCF \
   -O $FILTERED_INDELS_VCF \
   --filter-name "indel_QD_lt_5" \
   --filter-expression "QD < 5.0" \
   --filter-name "snp_FS_gt_55" \
   --filter-expression "FS > 55.0" \
   --filter-name "snp_ReadPosRankSum_lt_-20" \
   --filter-expression "ReadPosRankSum < -20.0" \
   --filter-name "snp_SQR_gt_10" \
   --filter-expression "SOR > 10.0" &

wait

#merge back together
ls cohort_filtered_*.vcf >merge.list

$SINGULARITY gatk MergeVcfs \
    -I merge.list \
    -O $MERGED_VARIANTS_VCF 


#$SINGULARITY gatk SelectVariants \
#     -R $REFERENCE \
#     -V cohort_filtered_SNPs_r3.vcf  \
#     -L AMPZ01026399.1 \
#     -O mito_variants.vcf

#create list of individuals with a lot of missing data
vcftools --remove-filtered all --min-alleles 2 --max-alleles 2 --max-missing 0.51 --mac 3 --vcf cohort_filtered_SNPs_r3_env.vcf --out cohort_filtered_SNPS_biallelic --recode --recode-INFO-all


#create list of individuals with a lot of missing data
vcftools --vcf cohort_filtered_SNPS_biallelic --missing-indv

cat out.imiss | awk '$5 >0.5 (print $1}' >data_poor_indivs.list

vcftools --vcf cohort_filtered_SNPS_biallelic.vcf --remove data_poor_indivs.list --recode --recode-INFO-all --out cohort_filtered_SNPS_biallelic_gt50P.vcf

#high proportion called variants:
vcftools --vcf cohort_filtered_SNPS_biallelic_gt50P.vcf.vcf --max-missing 0.95 --maf 0.05 --recode --recode-INFO-all --out cohort_filtered_SNPS_biallelic_gt50P_maf05_miss95 --min-meanDP 10

--thin


--SNPdensity <integer>




