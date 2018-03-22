#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr_recalibration.sh - use hard filtered SNPs to recalibrate base calls (round 1)

# TODO(nplatt): update comments

cd $BQSR_DIR

#!!! UNMODIFIED_BAM_DIR set in the pipeline script
for BAM in $(ls $UNMODIFIED_BAM_DIR/*.bam); do

    SAMPLE=$(echo $(basename $BAM) | cut -f1,2 -d"_")

    # BUILD TABLE FROM HQ SNPS AND PREV BAMS -----------------------------------    
    JOB_NAME="snp.$SAMPLE.recal_score_reads_"$ROUND
    THREADS=12
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid snp.merge_filtered_variants"
    SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    IN_BAM=$BAM
    IN_VCF=$BQSR_DIR/$ROUND"_vcfs/cohort_"$ROUND"_filteredVariants.g.vcf"
    OUT_TABLE=$BQSR_DIR/$ROUND"_tables/"$SAMPLE"_PRErecal."$ROUND".table"
    
    CMD="$SINGULARITY gatk BaseRecalibrator \
        -R $REFERENCE \
        -I $IN_BAM \
        --known-sites $IN_VCF \
        -O $OUT_TABLE"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"

    # MODIFY THE BAM FILES -----------------------------------------------------    
    JOB_NAME="snp.$SAMPLE.recal_modify_reads_$ROUND"
    THREADS=12
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid snp.$SAMPLE.recal_score_reads_$ROUND"
    SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

    IN_BAM=$IN_BAM
    IN_VCF=$IN_VCF
    IN_TABLE=$OUT_TABLE
    OUT_BAM=$BQSR_DIR/$ROUND"_bqsr_bams"/$SAMPLE"_bqsr-"$ROUND".bam"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"
    
    CMD="$SINGULARITY gatk ApplyBQSR \
        -R $REFERENCE \
        -I $IN_BAM \
        --bqsr-recal-file $IN_TABLE \
        -O $OUT_BAM"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB" 

    # SCORE THE RECALIBRATED BAM FILES -----------------------------------------    

    JOB_NAME="snp.$SAMPLE.recal_score_modreads_$ROUND"
    THREADS=12
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid snp.$SAMPLE.recal_modify_reads_$ROUND"
    SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

    IN_BAM=$OUT_BAM
    IN_VCF=$IN_VCF
    OUT_TABLE=$BQSR_DIR/$ROUND"_tables/"$SAMPLE"_POSTrecal.$ROUND.table"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"
    
    CMD="$SINGULARITY gatk BaseRecalibrator \
        -R $REFERENCE \
        -I $IN_BAM \
        --known-sites $IN_VCF \
        -O $OUT_TABLE"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB" 

    # ANALYZE THE RECALIBRATED DATA --------------------------------------------    
    JOB_NAME="snp.$SAMPLE.covariate_$ROUND"
    THREADS=12
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid snp.$SAMPLE.recal_score_reads_$ROUND,snp.$SAMPLE.recal_score_modreads_$ROUND"
    SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    #!!! BEFORE_MOD_TABLE set in the pipeline script
    BEFORE_MOD_TABLE=$PREV_RECAL_TABLE_DIR$SAMPLE$PREV_RECAL_TABLE_EXT
    AFTER_MOD_TABLE=$BQSR_DIR/$ROUND"_tables/"$SAMPLE"_POSTrecal.$ROUND.table"
    OUT_PDF=$BQSR_DIR/$ROUND"_tables/"$SAMPLE"_recalibration_plot."$ROUND".pdf"

    CMD="$SINGULARITY gatk AnalyzeCovariates \
        -before $BEFORE_MOD_TABLE \
        -after $AFTER_MOD_TABLE \
        -plots $OUT_PDF"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB" 

done


# ANALYZE ALL RECAL DATA -------------------------------------------------------    
PRE_LIST=all_PRE_$ROUND.list
POST_LIST=all_POST_$ROUND.list
PRE_TABLE=all_PRE_$ROUND.table
POST_TABLE=all_POST_$ROUND.table
ALL_PDF=all_recalibration_plot.$ROUND.pdf

ls *"_PRErecal.$ROUND.table" >$PRE_LIST
ls *"_POSTrecal.$ROUND.table" >$POST_LIST

$SINGULARITY gatk GatherBQSRReports --input $PRE_LIST --output $PRE_TABLE &
$SINGULARITY gatk GatherBQSRReports --input $POST_LIST --output $POST_TABLE &
wait
$SINGULARITY gatk AnalyzeCovariates -before $PRE_TABLE -after $POST_TABLE -plots $ALL_PDF

#check for convergence in recal data - analysis done directly

