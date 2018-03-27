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
#2) sites that are NOT snps or indels


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



#create list of individuals with a lot of missing data
vcftools \
    --remove-filtered-all \
    --min-alleles 2 \
    --max-alleles 2 \
    --max-missing 0.51 \
    --mac 3 \
    --vcf cohort_filtered_variants_r3.vcf \
    --out cohort_filtered_variants_biallelic \
    --recode \
    --recode-INFO-all


#create list of individuals with a lot of missing data
vcftools \
    --vcf cohort_filtered_variants_biallelic.recode.vcf \
    --missing-indv \

#remove individuals missing more than 25% of the data    
cat out.imiss | awk '$5 >0.25 {print $1}' >data_poor_indivs.list

vcftools --vcf cohort_filtered_variants_biallelic.recode.vcf --remove data_poor_indivs.list --recode --recode-INFO-all --out cohort_filtered_SNPS_biallelic_gt75P

#remove sites where 95% are genotyped
#high proportion called variants:
vcftools \
    --vcf cohort_filtered_SNPS_biallelic_gt75P.recode.vcf \
    --max-missing 0.95 \
    --maf 0.025 \
    --recode \
    --recode-INFO-all \
    --out cohort_filtered_SNPS_biallelic_gt75P_miss95 \
    --min-meanDP 10


#remove snps close to indels and clumps of indels
bcftools filter \
    --SnpGap 50 \
    --IndelGap 100 \
    --output cohort_filtered_SNPS_biallelic_gt75P_miss95_100bpfilt.vcf \
    -O v \
    cohort_filtered_SNPS_biallelic_gt75P_miss95.recode.vcf

#add unique identifer to ID field (replace "."
bcftools annotate --set-id +'%CHROM\:%POS' cohort_filtered_SNPS_biallelic_gt75P_miss95_100bpfilt.vcf >sHaem_filtered.vcf




