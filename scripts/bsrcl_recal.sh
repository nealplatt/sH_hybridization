#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bsrcl_recal.sh - use hard filtered SNPs to recalibrate base calls (round 1)

# TODO(nplatt): update comments

source master/nplatt/sH_hybridizationscripts/set_env.sh

cd $BSRCL_DIR

# SELECT SNPS ------------------------------------------------------------------
SELECT_SNPS_JOB_NAME=cohort_select_snps
THREADS=1

IN_VCF="$BSRCL_DIR/cohort_raw.g.vcf"
OUT_VCF="$BSRCL_DIR/cohort_raw_snps.g.vcf"
    
SELECT_SNPS="$SINGULARITY gatk SelectVariants -V $IN_VCF -select-type SNP -O $OUT_VCF -R $REFERENCE"

SELECT_SNPS_QSUB="$QSUB -pe mpi $THREADS -N $SELECT_SNPS_JOB_NAME -o logs/$SELECT_SNPS_JOB_NAME.log"

echo $SELECT_SNPS >scripts/$SELECT_SNPS_QSUB.sh

cat scripts/$SELECT_SNPS_JOB_NAME.sh | $SELECT_SNPS_QSUB

# FILTER SNPS ------------------------------------------------------------------
FILTER_SNPS_JOB_NAME=cohort_FILTER_snps
THREADS=1

IN_VCF=$OUT_VCF
OUT_VCF="$BSRCL_DIR/cohort_filtered_snps.g.vcf"
    
FILTER_SNPS="$SINGULARITY gatk VariantFiltration \
    -R $REFERENCE \
    -V $IN_VCF \
    --filterExpression "'"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"'" \
    --filterName "'"recommended_snp_filter"'" \
    -o $OUT_VCF"

FILTER_SNPS_QSUB="$QSUB -pe mpi $THREADS -N $FILTER_SNPS_JOB_NAME -o logs/$FILTER_SNPS_JOB_NAME.log -hold_jid $SELECT_SNPS_JOB_NAME"

echo $FILTER_SNPS >scripts/$FILTER_SNPS_QSUB.sh

cat scripts/$FILTER_SNPS_JOB_NAME.sh | $FILTER_SNPS_QSUB

# SELECT INDELS ----------------------------------------------------------------
SELECT_INDELS_JOB_NAME=cohort_select_indels
THREADS=1

IN_VCF="$BSRCL_DIR/cohort_raw.g.vcf"
OUT_VCF="$BSRCL_DIR/cohort_raw_indels.g.vcf"
    
SELECT_INDELS="$SINGULARITY gatk SelectVariants -V $IN_VCF -select-type INDEL -O $OUT_VCF -R $REFERENCE"

SELECT_INDELS_QSUB="$QSUB -pe mpi $THREADS -N $SELECT_INDELS_JOB_NAME -o logs/$SELECT_INDELS_JOB_NAME.log"

echo $SELECT_INDELS >scripts/$SELECT_INDELS_QSUB.sh

cat scripts/$SELECT_JOB_NAME.sh | $SELECT_INDELS_QSUB


# FILTER INDELS ------------------------------------------------------------------
FILTER_INDELS_JOB_NAME=cohort_FILTER_INDELS
THREADS=1

IN_VCF=$OUT_VCF
OUT_VCF="$BSRCL_DIR/cohort_filtered_indels.g.vcf"
    
FILTER_INDELS="$SINGULARITY gatk VariantFiltration \
    -R $REFERENCE \
    -V $IN_VCF \
    --filterExpression "'"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"'" \
    --filterName "'"recommended_indel_filter"'" \
    -o $OUT_VCF"

FILTER_INDELS_QSUB="$QSUB -pe mpi $THREADS -N $FILTER_INDELS_JOB_NAME -o logs/$FILTER_INDELS_JOB_NAME.log -hold_jid $SELECT_INDELS_JOB_NAME"

echo $FILTER_INDELS >scripts/$FILTER_INDELS_QSUB.sh

cat scripts/$FILTER_INDELS_JOB_NAME.sh | $FILTER_INDELS_QSUB


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
