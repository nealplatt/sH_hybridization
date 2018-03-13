#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr-r2_bqsr.sh - use hard filtered SNPs to recalibrate base calls (round 1)

# TODO(nplatt): update comments

source /master/nplatt/sH_hybridization/scripts/set-env.sh

mkdir $BQSR_DIR/bqsr_r2 $BQSR_DIR/r2_cov_plots $BQSR_DIR/r2_bqsr_bams

cd $BQSR_DIR/bqsr_r2

# SELECT SNPS ------------------------------------------------------------------
JOB_NAME=cohort_select_snps
THREADS=1
LOG="$LOGS_DIR/$JOB_NAME.log" 
DEPEND=""
SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

IN_VCF="$BQSR_DIR/cohort_r1_BQSR.g.vcf"
OUT_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_BSQR_rawSNPS.g.vcf"
    
CMD="$SINGULARITY gatk SelectVariants \
    -V $IN_VCF \
    -select-type SNP \
    -O $OUT_VCF \
    -R $REFERENCE"

DELETE $LOG $SCRIPT
SUBMIT $CMD $SCRIPT $JOB_QSUB

# FILTER SNPS ------------------------------------------------------------------
JOB_NAME="cohort_filter_snps"
THREADS=1
LOG="$LOGS_DIR/$JOB_NAME.log" 
DEPEND="-hold_jid cohort_select_snps"
SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

IN_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_BSQR_rawSNPS.g.vcf"
OUT_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_filtered_snps.g.vcf"

CMD="$SINGULARITY gatk VariantFiltration \
    -R $REFERENCE \
    -V $IN_VCF \
    --filter-expression "'"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"'" \
    --filter-name "'"recommended_snp_filter"'" \
    -O $OUT_VCF"

DELETE $LOG $SCRIPT
SUBMIT $CMD $SCRIPT $JOB_QSUB

# SELECT INDELS ----------------------------------------------------------------
JOB_NAME="cohort_select_indels"
THREADS=1
LOG="$LOGS_DIR/$JOB_NAME.log" 
DEPEND=""
SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

IN_VCF="$BQSR_DIR/cohort_r1_BQSR.g.vcf"
OUT_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_BSQR_rawINDELS.g.vcf"
    
CMD="$SINGULARITY gatk SelectVariants \
    -V $IN_VCF \
    -select-type INDEL \
    -O $OUT_VCF \
    -R $REFERENCE"

DELETE $LOG $SCRIPT
SUBMIT $CMD $SCRIPT $JOB_QSUB


# FILTER INDELS ------------------------------------------------------------------
JOB_NAME="cohort_filter_indels"
THREADS=1
LOG="$LOGS_DIR/$JOB_NAME.log" 
DEPEND="-hold_jid cohort_select_indels"
SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

IN_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_BSQR_rawINDELS.g.vcf"
OUT_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_filtered_INDELS.g.vcf"

CMD="$SINGULARITY gatk VariantFiltration \
    -R $REFERENCE \
    -V $IN_VCF \
    --filter-expression "'"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"'" \
    --filter-name "'"recommended_indel_filter"'" \
    -O $OUT_VCF"

DELETE $LOG $SCRIPT
SUBMIT $CMD $SCRIPT $JOB_QSUB


################################################################################
#                       R E C A L I B R A T I O N 
################################################################################

for BAM in $(ls $MAP_DIR/*_processed.bam); do

    SAMPLE=$(basename $BAM _processed.bam)

    # RECAL SCORE/OBSERVE READS-------------------------------------------------    
    JOB_NAME="$SAMPLE.recal_score_reads_r2"
    THREADS=12
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid merge_filtered_vcfs"
    SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    IN_BAM=$BAM
    IN_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_filtered_snps.g.vcf"
    OUT_TABLE=$BQSR_DIR/r2_cov_plots/$SAMPLE"_recal-1_data.table"
    
    CMD="$SINGULARITY gatk BaseRecalibrator \
        -R $REFERENCE \
        -I $IN_BAM \
        --known-sites $IN_VCF \
        -O $OUT_TABLE"

    #DELETE $LOG $SCRIPT
    #SUBMIT $CMD $SCRIPT $JOB_QSUB

    # MOD READS ----------------------------------------------------------------    
    JOB_NAME="$SAMPLE.recal_modify_reads_r2"
    THREADS=12
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid $SAMPLE.recal_score_reads_r2"
    SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

    IN_BAM=$BAM
    IN_VCF=$IN_VCF
    IN_TABLE=$OUT_TABLE
    OUT_BAM=$BQSR_DIR/r2_bqsr_bams/$SAMPLE".bqsr-2.bam"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"
    
    CMD="$SINGULARITY gatk ApplyBQSR \
        -R $REFERENCE \
        -I $IN_BAM \
        --bqsr-recal-file $IN_TABLE \
        -O $OUT_BAM"

    #DELETE $LOG $SCRIPT
    #SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB" 

    # RECAL SCORE/OBSERVE READS-------------------------------------------------    

    JOB_NAME="$SAMPLE.recal_score_modreads_r2"
    THREADS=12
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid $SAMPLE.recal_modify_reads_r2"
    SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

    IN_BAM=$OUT_BAM
    IN_VCF=$IN_VCF
    OUT_TABLE=$BQSR_DIR/r2_cov_plots/$SAMPLE"_postrecal-1_data.table"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"
    
    CMD="$SINGULARITY gatk BaseRecalibrator \
        -R $REFERENCE \
        -I $IN_BAM \
        --known-sites $IN_VCF \
        -O $OUT_TABLE"

    #DELETE $LOG $SCRIPT
    #SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB" 

    # RECAL ANALYZE COVARIATES -------------------------------------------------    
    JOB_NAME="$SAMPLE.covariate_r2"
    THREADS=12
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid $SAMPLE.recal_score_reads_r2,$SAMPLE.recal_score_modreads_r2"
    SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    BEFORE_TABLE=$BQSR_DIR/r1_cov_plots/$SAMPLE"_postrecal-1_data.table"
    AFTER_TABLE=$BQSR_DIR/r2_cov_plots/$SAMPLE"_postrecal-1_data.table"
    OUT_PDF=$BQSR_DIR/r2_cov_plots/$SAMPLE"_recalibration_plot.2.pdf"

    CMD="$SINGULARITY gatk AnalyzeCovariates \
        -before $BEFORE_TABLE \
        -after $AFTER_TABLE \
        -plots $OUT_PDF"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB" 

done

#check for convergence in recal data

# CLEANUP ----------------------------------------------------------------------
rm -r hc_vcf_r2
rm -r individual_vcf_r2
rm -r r1_bqsr_bams
rm -r db_r2
rm -r genotype_interval_vcfs_r2
rm -r cohort_vcf_r2
rm samples_r2.list
rm samples.list
rm cohort_filtered_snps.g.vcf.idx
rm cohort_preBSQR.g.vcf
rm cohort_preBSQR.g.vcf.idx
rm cohort_preBSQR_rawIndels.g.vcf
rm cohort_preBSQR_rawIndels.g.vcf.idx
rm cohort_preBSQR_rawSNPS.g.vcf
rm cohort_preBSQR_rawSNPS.g.vcf.idx
rm cohort_filtered_indels.g.vcf
rm cohort_filtered_indels.g.vcf.idx
rm cohort_filtered_snps.g.vcf
rm cohort_filtered_snps.g.vcf.idx
rm cohort_r1_BQSR.g.vcf



