#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bsrcl_recal.sh - use hard filtered SNPs to recalibrate base calls (round 1)

# TODO(nplatt): update comments

source master/nplatt/sH_hybridizationscripts/set_env.sh


#select all of the SNPs
singularity exec ~/snpCalling_v0.0.5.img gatk SelectVariants -V cohort_sorted.g.vcf -select-type SNP -O cohort_raw_snps.g.vcf -R /master/nplatt/sH_hybridization/data/genome/schHae_v1.fa

#select all of the indels
singularity exec ~/snpCalling_v0.0.5.img gatk SelectVariants -V cohort_sorted.g.vcf -select-type INDEL -O cohort_raw_indels.g.vcf -R /master/nplatt/sH_hybridization/data/genome/schHae_v1.fa

#filter out the low qual SNPs
singularity exec ~/snpCalling_v0.0.5.img gatk VariantFiltration \
    -R -R /master/nplatt/sH_hybridization/data/genome/schHae_v1.fa \
    -V cohort_raw_snps.g.vcf \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "initial_snp_filter" \
    -o filtered_snps.vcf 


################################################################################
#recalibration
################################################################################
#1

#(add read groups)
#(sort)
#(mark duplicates)


REF=
for BAM in LIST_OF_SAMPLES; do
    RECAL_1=singularity exec snpCalling_v0.0.5.img gatk BaseRecalibrator \
        -R /master/nplatt/sH_hybridization/data/genome/schHae_v1.fa \
        -I $BAM_FILE \
        -knownSites filtered_snps.vcf \
        -o <sample>_recal_data.1.table

    RECAL_2=singularity exec snpCalling_v0.0.5.img gatk BaseRecalibrator \
        -R $REF \
        -I $BAM_FILE \
        -knownSites recal_0.filtered.g.vcf \
        -BQSR recal_data.1.table \
        -o post_recal_data.1.table

    ANA_COV=singularity exec snpCalling_v0.0.5.img gatk AnalyzeCovariates \
        -R $REF \
        -before recal_data.1.table \
        -after post_recal_data.1.table \
        -plots recalibration_plots.1.pdf


    submit recal-1
    submit recal-2 -hold recal1
    submit ana_cov -hold recal2


done
#2

1.1.4 HaplotypeCaller -R $REF -I $BAM_FILE -ERC GVCF -BQSR post_recal_data.1.table -o recal_1.g.vcf
1.1.5 Filter recal_1.g.vcf --> recal_1.g.filtered.vcf

Iteration 2
1.2.1 BaseRecalibrator -R $REF -I $BAM_FILE -knownSites recal_1.g.filtered.vcf -o recal_data.2.table
1.2.2 BaseRecalibrator -R $REF -I $BAM_FILE -knownSites recal_1.g.filtered.vcf -BQSR recal_data.2.table -o post_recal_data.2.table
1.2.3 AnalyzeCovariates -R $REF -before recal_data.2.table -after post_recal_data.2.table -plots recalibration_plots.2.pdf
1.2.4 HaplotypeCaller -R $REF -I $BAM_FILE -ERC GVCF -BQSR post_recal_data.2.table -o recal_2.g.vcf
1.2.5 Filter recal_2.g.vcf --> recal_2.g.filtered.vcf

--> HC -BQSR on the fly recalibration of "raw" bam files
--> Each iteration should produce a new vcf file that is used for the recalibration of "raw" reads
--> Differences between the AnalyzeCovariates before and after should diminush after each iteration (difference 1.1.3 > difference 1.2.3)
