#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bsrcl_recal.sh - use hard filtered SNPs to recalibrate base calls (round 1)

# TODO(nplatt): update comments
# TODO(nplatt): has not been through test-run yet

source master/nplatt/sH_hybridizationscripts/set_env.sh

cd $BSRCL_DIR

# SELECT SNPS ------------------------------------------------------------------
SELECT_SNPS_JOB_NAME=cohort_select_snps
THREADS=1

IN_VCF="$BSRCL_DIR/cohort_preBSQR.g.vcf"
OUT_VCF="$BSRCL_DIR/cohort_preBSQR_rawSNPS.g.vcf"
    
SELECT_SNPS="$SINGULARITY gatk SelectVariants \
    -V $IN_VCF \
    -select-type SNP \
    -O $OUT_VCF \
    -R $REFERENCE"

SELECT_SNPS_QSUB="$QSUB -pe mpi $THREADS -N $SELECT_SNPS_JOB_NAME -o logs/$SELECT_SNPS_JOB_NAME.log"

echo $SELECT_SNPS >scripts/$SELECT_SNPS_JOB_NAME.sh

cat scripts/$SELECT_SNPS_JOB_NAME.sh | $SELECT_SNPS_QSUB

# FILTER SNPS ------------------------------------------------------------------
FILTER_SNPS_JOB_NAME="cohort_filter_snps"
THREADS=1

IN_VCF="$BSRCL_DIR/cohort_preBSQR_rawSNPS.g.vcf"
OUT_VCF="$BSRCL_DIR/cohort_filtered_snps.g.vcf"

FILTER_SNPS="$SINGULARITY gatk VariantFiltration \
    -R $REFERENCE \
    -V $IN_VCF \
    --filter-expression "'"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"'" \
    --filter-name "'"recommended_snp_filter"'" \
    -O $OUT_VCF"

FILTER_SNPS_QSUB="$QSUB -pe mpi $THREADS -N $FILTER_SNPS_JOB_NAME -o logs/$FILTER_SNPS_JOB_NAME.log -hold_jid $SELECT_SNPS_JOB_NAME"

echo $FILTER_SNPS >scripts/$FILTER_SNPS_JOB_NAME.sh

cat scripts/$FILTER_SNPS_JOB_NAME.sh | $FILTER_SNPS_QSUB

# SELECT INDELS ----------------------------------------------------------------
SELECT_INDELS_JOB_NAME="cohort_select_indels"
THREADS=1

IN_VCF="$BSRCL_DIR/cohort_preBSQR.g.vcf"
OUT_VCF="$BSRCL_DIR/cohort_preBSQR_rawIndels.g.vcf"
    
SELECT_INDELS="$SINGULARITY gatk SelectVariants \
    -V $IN_VCF \
    -select-type INDEL \
    -O $OUT_VCF \
    -R $REFERENCE"

SELECT_INDELS_QSUB="$QSUB -pe mpi $THREADS -N $SELECT_INDELS_JOB_NAME -o logs/$SELECT_INDELS_JOB_NAME.log"

echo $SELECT_INDELS >scripts/$SELECT_JOB_NAME.sh

cat scripts/$SELECT_JOB_NAME.sh | $SELECT_INDELS_QSUB


# FILTER INDELS ------------------------------------------------------------------
FILTER_INDELS_JOB_NAME="cohort_filter_indels"
THREADS=1

IN_VCF="$BSRCL_DIR/cohort_preBSQR_rawIndels.g.vcf"
OUT_VCF="$BSRCL_DIR/cohort_filtered_indels.g.vcf"

FILTER_INDELS="$SINGULARITY gatk VariantFiltration \
    -R $REFERENCE \
    -V $IN_VCF \
    --filter-expression "'"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"'" \
    --filter-name "'"recommended_indel_filter"'" \
    -O $OUT_VCF"

FILTER_INDELS_QSUB="$QSUB -pe mpi $THREADS -N $FILTER_INDELS_JOB_NAME -o logs/$FILTER_INDELS_JOB_NAME.log -hold_jid $SELECT_INDELS_JOB_NAME"

echo $FILTER_INDELS >scripts/$FILTER_INDELS_JOB_NAME.sh

cat scripts/$FILTER_INDELS_JOB_NAME.sh | $FILTER_INDELS_QSUB


# MERGE FILTERED VCFS ------------------------------------------------------------------
MERGE_JOB_NAME="merge_filtered_vcfs"
THREADS=1

IN_SNP_VCF="$BSRCL_DIR/cohort_filtered_indels.g.vcf"
IN_INDEL_VCF="$BSRCL_DIR/cohort_filtered_snps.g.vcf"
OUT_VCF="$BSRCL_DIR/cohort_filtered_snps-indels.g.vcf"
        
MERGE="$SINGULARITY gatk CombineVariants \
   -R $REFERENCE \
   --variant:snps $IN_SNP_VCF \
   --variant:indels $IN_INDEL_VCF \
   -O $OUT_VCF \
   -genotypeMergeOptions PRIORITIZE \
   -priority snps,indels"

MERGE_QSUB="$QSUB -pe mpi $THREADS -N $MERGE_JOB_NAME -o logs/$MERGE_JOB_NAME.log -hold_jid $FILTER_INDELS_JOB_NAME,$FILTER_SNPS_JOB_NAME"

echo $MERGE >scripts/$MERGE_JOB_NAME.sh

cat scripts/$MERGE_JOB_NAME.sh | $MERGE_QSUB


################################################################################
#                       R E C A L I B R A T I O N 
################################################################################

for BAM in $(ls $MAP_DIR/*_processed.bam); do

    SAMPLE=$(basename $BAM _processed.bam)

    # RECAL SCORE/OBSERVE READS-------------------------------------------------    
    RECAL_OBSERVE_JOB_NAME="$SAMPLE.recal_score_reads"
    THREADS=12

    IN_BAM=$MAP_DIR/$BAM
    IN_VCF="$BSRCL_DIR/cohort_filtered_snps.g.vcf"
    OUT_TABLE=$SAMPLE"_recal-1_data.table"
    
    RECAL_OBSERVE="$SINGULARITY gatk BaseRecalibrator \
        -R $REFERENCE \
        -I $BAM \
        --known-sites $IN_VCF \
        -O $OUT_TABLE"

    RECAL_OBSERVE_QSUB="$QSUB -pe mpi $THREADS -N $RECAL_OBSERVE_JOB_NAME -o logs/$RECAL_OBSERVE_JOB_NAME.log -hold_jid $FILTER_SNPS_JOB_NAME"

    #echo $RECAL_OBSERVE >scripts/$RECAL_OBSERVE_JOB_NAME.sh
    #cat scripts/$RECAL_OBSERVE_JOB_NAME.sh | $RECAL_OBSERVE_QSUB

    # MOD READS ----------------------------------------------------------------    
    MOD_READS_JOB_NAME="$SAMPLE.recal_modifed_reads"
    THREADS=12

    IN_BAM=$BAM
    IN_VCF=$IN_VCF
    IN_TABLE=$OUT_TABLE


    OUT_BAM=$SAMPLE".bqsr-1.bam"
    
    RECAL_MODIFY="$SINGULARITY gatk ApplyBQSR \
        -R $REFERENCE \
        -I $IN_BAM \
        --bqsr-recal-file $IN_TABLE \
        -O $OUT_BAM"

    RECAL_MODIFY_QSUB="$QSUB -pe mpi $THREADS -N $RECAL_MODIFY_JOB_NAME -o logs/$RECAL_MODIFY_JOB_NAME.log, -hold_jid $RECAL_OBSERVE_JOB_NAME"

    #echo $RECAL_MODIFY >scripts/$RECAL_MODIFY_JOB_NAME.sh
    #cat scripts/$RECAL_MODIFY_JOB_NAME.sh | $RECAL_MODIFY_QSUB    

    # RECAL SCORE/OBSERVE READS-------------------------------------------------    

    SCORE_MOD_JOB_NAME="$SAMPLE.recal_score_reads"
    THREADS=12

    IN_BAM=$OUT_BAM
    IN_VCF="$BSRCL_DIR/cohort_filtered_snps.g.vcf"
    OUT_TABLE=$SAMPLE"_postrecal-1_data.table"
    
    SCORE_MOD="$SINGULARITY gatk BaseRecalibrator \
        -R $REFERENCE \
        -I $BAM \
        --known-sites $IN_VCF \
        -O $OUT_TABLE"

    SCORE_MOD_QSUB="$QSUB -pe mpi $THREADS -N $SCORE_MOD_JOB_NAME -o logs/$SCORE_MOD_JOB_NAME.log -hold_jid $MOD_READS_JOB_NAME"

    #echo $SCORE_MOD >scripts/$SCORE_MOD_JOB_NAME.sh
    #cat scripts/$SCORE_MOD_JOB_NAME.sh | $SCORE_MOD_QSUB

    # RECAL ANALYZE COVARIATES -------------------------------------------------    
    ANALYZE_JOB_NAME="$SAMPLE.covariate_analysis"
    THREADS=12

    BEFORE_TABLE=$SAMPLE"_recal-1_data.table"
    AFTER_TABLE=$SAMPLE"_postrecal-1_data.table"
    OUT_PDF=$SAMPLE"_recalibration_plot.1.pdf"

    ANALYZE="$SINGULARITY gatk AnalyzeCovariates \
        -before $BEFORE_TABLE \
        -after $AFTER_TABLE \
        -plots $OUT_PDF"

    ANALYZE_QSUB="$QSUB -pe mpi $THREADS -N $ANALYZE_JOB_NAME -o logs/$ANALYZE_JOB_NAME.log -hold_jid $SCORE_MOD_JOB_NAME,$RECAL_OBSERVE_JOB_NAME"

    echo $ANALYZE >scripts/$ANALYZE_JOB_NAME.sh
    cat scripts/$ANALYZE_JOB_NAME.sh | $ANALYZE_QSUB  

done
