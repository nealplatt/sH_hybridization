#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr-r2_bqsr.sh - use hard filtered SNPs to recalibrate base calls (round 1)

# TODO(nplatt): update comments

cd $BQSR_DIR


for BAM in $(ls $MAP_DIR/*_processed.bam); do <-----------------------------------------------------changed to in.bam from prev round

    SAMPLE=$(basename $BAM _processed.bam)<-----------------------------------------------------changed to in.bam from prev round

    # RECAL SCORE/OBSERVE READS-------------------------------------------------    
    JOB_NAME="snp.$SAMPLE.recal_score_reads_"$ROUND
    THREADS=12
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid snp.merge_filtered_variants"
    SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    IN_BAM=$BAM
    IN_VCF=$ROUND"_vcfs/cohort_"$ROUND"_filteredVariants.g.vcf"
    OUT_TABLE=$BQSR_DIR/$ROUND"_tables/"$SAMPLE"_prerecal-"$ROUND".table"
    
    CMD="$SINGULARITY gatk BaseRecalibrator \
        -R $REFERENCE \
        -I $IN_BAM \
        --known-sites $IN_VCF \
        -O $OUT_TABLE"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"

    # MOD READS ----------------------------------------------------------------    
    JOB_NAME="snp.$SAMPLE.recal_modify_reads_r2"
    THREADS=12
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid snp.$SAMPLE.recal_score_reads_r2"
    SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

    IN_BAM=$BAM
    IN_VCF=$IN_VCF
    IN_TABLE=$OUT_TABLE
    OUT_BAM=$BQSR_DIR/$ROUND"_bqsr_bams"/$SAMPLE".bqsr-"$ROUND".bam"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"
    
    CMD="$SINGULARITY gatk ApplyBQSR \
        -R $REFERENCE \
        -I $IN_BAM \
        --bqsr-recal-file $IN_TABLE \
        -O $OUT_BAM"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB" 

    # RECAL SCORE/OBSERVE READS-------------------------------------------------    

    JOB_NAME="snp.$SAMPLE.recal_score_modreads_r2"
    THREADS=12
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid snp.$SAMPLE.recal_modify_reads_r2"
    SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

    IN_BAM=$OUT_BAM
    IN_VCF=$IN_VCF
---------------------------------------------OUT_TABLE=$BQSR_DIR/r2_cov_plots/$SAMPLE"_postrecal-1_data.table"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"
    
    CMD="$SINGULARITY gatk BaseRecalibrator \
        -R $REFERENCE \
        -I $IN_BAM \
        --known-sites $IN_VCF \
        -O $OUT_TABLE"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB" 

    # RECAL ANALYZE COVARIATES -------------------------------------------------    
    JOB_NAME="snp.$SAMPLE.covariate_r2"
    THREADS=12
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid snp.$SAMPLE.recal_score_reads_r2,snp.$SAMPLE.recal_score_modreads_r2"
    SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

-----------------------------------------BEFORE_TABLE=$BQSR_DIR/r1_cov_plots/$SAMPLE"_postrecal-1_data.table"
-----------------------------------------AFTER_TABLE=$BQSR_DIR/r2_cov_plots/$SAMPLE"_postrecal-1_data.table"
    OUT_PDF=$BQSR_DIR/$ROUND"_cov_plots/"$SAMPLE"_recalibration_plot."$ROUND".pdf"

    CMD="$SINGULARITY gatk AnalyzeCovariates \
        -before $BEFORE_TABLE \
        -after $AFTER_TABLE \
        -plots $OUT_PDF"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB" 

done

#check for convergence in recal data

# CLEANUP ----------------------------------------------------------------------
#rm -r hc_vcf_r2
#rm -r individual_vcf_r2
#rm -r r1_bqsr_bams
#rm -r db_r2
#rm -r genotype_interval_vcfs_r2
#rm -r cohort_vcf_r2
#rm samples_r2.list
#rm samples.list
#rm cohort_filtered_snps.g.vcf.idx
#rm cohort_preBSQR.g.vcf
#rm cohort_preBSQR.g.vcf.idx
#rm cohort_preBSQR_rawIndels.g.vcf
#rm cohort_preBSQR_rawIndels.g.vcf.idx
#rm cohort_preBSQR_rawSNPS.g.vcf
#rm cohort_preBSQR_rawSNPS.g.vcf.idx
#rm cohort_filtered_indels.g.vcf
#rm cohort_filtered_indels.g.vcf.idx
#rm cohort_filtered_snps.g.vcf
#rm cohort_filtered_snps.g.vcf.idx
#rm cohort_r1_BQSR.g.vcf



