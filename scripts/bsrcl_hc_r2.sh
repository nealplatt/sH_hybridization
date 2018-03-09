#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bsrcl_hc_r2.sh - use BSQR tables for second round of HC in recalibration proc.

# TODO(nplatt): update comments
# TODO(nplatt): has not been through test-run yet

source master/nplatt/sH_hybridizationscripts/set_env.sh

cd $BSRCL_DIR

#-------------------------------------------------------------------------------
#NOTES FROM GATK FORUM:
#https://gatkforums.broadinstitute.org/gatk/discussion/4913/recommended-protocol-for-bootstrapping-haplotypecaller-and-baserecalibrator-outputs

#1.1.4 HaplotypeCaller -R $REF -I $BAM_FILE -ERC GVCF -BQSR post_recal_data.1.table -o recal_1.g.vcf
#1.1.5 Filter recal_1.g.vcf --> recal_1.g.filtered.vcf

#Iteration 2
#1.2.1 BaseRecalibrator -R $REF -I $BAM_FILE -knownSites recal_1.g.filtered.vcf -o recal_data.2.table
#1.2.2 BaseRecalibrator -R $REF -I $BAM_FILE -knownSites recal_1.g.filtered.vcf -BQSR recal_data.2.table -o post_recal_data.2.table
#1.2.3 AnalyzeCovariates -R $REF -before recal_data.2.table -after post_recal_data.2.table -plots recalibration_plots.2.pdf
#1.2.4 HaplotypeCaller -R $REF -I $BAM_FILE -ERC GVCF -BQSR post_recal_data.2.table -o recal_2.g.vcf
#1.2.5 Filter recal_2.g.vcf --> recal_2.g.filtered.vcf

#--> HC -BQSR on the fly recalibration of "raw" bam files
#--> Each iteration should produce a new vcf file that is used for the recalibration of "raw" reads
#--> Differences between the AnalyzeCovariates before and after should diminush after each iteration (difference 1.1.3 > difference 1.2.3)
