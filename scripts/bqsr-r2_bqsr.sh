#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr-r2_bqsr.sh - use hard filtered SNPs to recalibrate base calls (round 1)

# TODO(nplatt): update comments

source /master/nplatt/sH_hybridization/scripts/set-env.sh

mkdir $BQSR_DIR/bqsr_r2

cd $BQSR_DIR/bqsr_r2

# SELECT SNPS ------------------------------------------------------------------
JOB_NAME=cohort_select_snps
THREADS=1
LOG="-o logs/$JOB_NAME.log" 
DEPEND=""
SCRIPT="scripts/$JOB_NAME.sh"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME $LOG $DEPEND"

IN_VCF="$BQSR_DIR/cohort_r1_BQSR.g.vcf"
OUT_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_BSQR_rawSNPS.g.vcf"
    
CMD="$SINGULARITY gatk SelectVariants \
    -V $IN_VCF \
    -select-type SNP \
    -O $OUT_VCF \
    -R $REFERENCE"

echo $CMD >$SCRIPT
cat $SCRIPT | $JOB_QSUB

# FILTER SNPS ------------------------------------------------------------------
JOB_NAME="cohort_filter_snps"
THREADS=1
LOG="-o logs/$JOB_NAME.log"
DEPEND="-hold_jid cohort_select_snps"
SCRIPT="scripts/$JOB_NAME.sh"

IN_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_BSQR_rawSNPS.g.vcf"
OUT_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_filtered_snps.g.vcf"

CMD="$SINGULARITY gatk VariantFiltration \
    -R $REFERENCE \
    -V $IN_VCF \
    --filter-expression "'"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"'" \
    --filter-name "'"recommended_snp_filter"'" \
    -O $OUT_VCF"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME $LOG $DEPEND"

echo $CMD >$SCRIPT
cat $SCRIPT | $JOB_QSUB

# SELECT INDELS ----------------------------------------------------------------
JOB_NAME="cohort_select_indels"
THREADS=1
LOG="-o logs/$JOB_NAME.log"
DEPEND=""
SCRIPT="scripts/$JOB_NAME.sh"

IN_VCF="$BQSR_DIR/cohort_r1_BQSR.g.vcf"
OUT_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_BSQR_rawINDELS.g.vcf"
    
CMD="$SINGULARITY gatk SelectVariants \
    -V $IN_VCF \
    -select-type INDEL \
    -O $OUT_VCF \
    -R $REFERENCE"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME $LOG $DEPEND"

echo $CMD >$SCRIPT
cat $SCRIPT | $JOB_QSUB


# FILTER INDELS ------------------------------------------------------------------
JOB_NAME="cohort_filter_indels"
THREADS=1
LOG="-o logs/$JOB_NAME.log"
DEPEND="-hold_jid cohort_select_indels"
SCRIPT="scripts/$JOB_NAME.sh"

IN_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_BSQR_rawINDELS.g.vcf"
OUT_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_filtered_INDELS.g.vcf"

CMD="$SINGULARITY gatk VariantFiltration \
    -R $REFERENCE \
    -V $IN_VCF \
    --filter-expression "'"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"'" \
    --filter-name "'"recommended_indel_filter"'" \
    -O $OUT_VCF"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME $LOG $DEPEND"

echo $CMD >$SCRIPT
cat $SCRIPT | $JOB_QSUB


# MERGE FILTERED VCFS ------------------------------------------------------------------
JOB_NAME="merge_filtered_vcfs"
THREADS=1
LOG="-o logs/$JOB_NAME.log"
DEPEND="-hold_jid cohort_filter_snps,cohort_filter_indels"
SCRIPT="scripts/$JOB_NAME.sh"

IN_SNP_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_filtered_SNPS.g.vcf"
IN_INDEL_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_filtered_INDELS.g.vcf"
OUT_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_filtered_variants.vcf"
        
CMD="$SINGULARITY gatk CombineVariants \
   -R $REFERENCE \
   --variant:snps $IN_SNP_VCF \
   --variant:indels $IN_INDEL_VCF \
   -O $OUT_VCF \
   -genotypeMergeOptions PRIORITIZE \
   -priority snps,indels"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME $LOG $DEPEND"

echo $CMD >$SCRIPT
cat $SCRIPT | $JOB_QSUB


################################################################################
#                       R E C A L I B R A T I O N 
################################################################################

for BAM in $(ls $MAP_DIR/*_processed.bam); do

    SAMPLE=$(basename $BAM _processed.bam)

    # RECAL SCORE/OBSERVE READS-------------------------------------------------    
    JOB_NAME="$SAMPLE.recal_score_reads_r2"
    THREADS=12
    LOG="-o logs/$JOB_NAME.log"
    DEPEND="-hold_jid merge_filtered_vcfs"
    SCRIPT="scripts/$JOB_NAME.sh"

    IN_BAM=$MAP_DIR/$BAM
    IN_VCF="$BQSR_DIR/bqsr_r2/cohort_r2_filtered_variants.vcf"
    OUT_TABLE=$BQSR_DIR/bqsr_r2/$SAMPLE"_recal-1_data.table"
    
    CMD="$SINGULARITY gatk BaseRecalibrator \
        -R $REFERENCE \
        -I $IN_BAM \
        --known-sites $IN_VCF \
        -O $OUT_TABLE"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME $LOG $DEPEND"

    echo $CMD >$SCRIPT
    cat $SCRIPT | $JOB_QSUB

    # MOD READS ----------------------------------------------------------------    
    JOB_NAME="$SAMPLE.recal_modifed_reads_r2"
    THREADS=12
    LOG="-o logs/$JOB_NAME.log"
    DEPEND="-hold_jid $SAMPLE.recal_score_reads"
    SCRIPT="scripts/$JOB_NAME.sh"

    IN_BAM=$BAM
    IN_VCF=$IN_VCF
    IN_TABLE=$OUT_TABLE
    OUT_BAM=$BQSR_DIR/bqsr_r2/$SAMPLE".bqsr-2.bam"
    
    CMD="$SINGULARITY gatk ApplyBQSR \
        -R $REFERENCE \
        -I $IN_BAM \
        --bqsr-recal-file $IN_TABLE \
        -O $OUT_BAM"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME $LOG $DEPEND"

    echo $CMD >$SCRIPT
    cat $SCRIPT | $JOB_QSUB  

    # RECAL SCORE/OBSERVE READS-------------------------------------------------    

    JOB_NAME="$SAMPLE.recal_score_modreads_r2"
    THREADS=12
    LOG="-o logs/$JOB_NAME.log"
    DEPEND="-hold_jid $SAMPLE.recal_score_modreads_r2"
    SCRIPT="scripts/$JOB_NAME.sh"

    IN_BAM=$OUT_BAM
    IN_VCF=$IN_VCF
    OUT_TABLE=$BQSR_DIR/bqsr_r2/$SAMPLE"_postrecal-1_data.table"
    
    CMD="$SINGULARITY gatk BaseRecalibrator \
        -R $REFERENCE \
        -I $IN_BAM \
        --known-sites $IN_VCF \
        -O $OUT_TABLE"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME $LOG $DEPEND"

    echo $CMD >$SCRIPT
    cat $SCRIPT | $JOB_QSUB 

    # RECAL ANALYZE COVARIATES -------------------------------------------------    
    JOB_NAME="$SAMPLE.covariate_r2"
    THREADS=12
    LOG="-o logs/$JOB_NAME.log"
    DEPEND="-hold_jid $SAMPLE.recal_score_reads_r2,$SAMPLE.recal_score_modreads_r2"
    SCRIPT="scripts/$JOB_NAME.sh"

    BEFORE_TABLE=$BQSR_DIR/bqsr_r2/$SAMPLE"_recal-1_data.table"
    AFTER_TABLE=$BQSR_DIR/bqsr_r2/$SAMPLE"_postrecal-1_data.table"
    OUT_PDF=$BQSR_DIR/bqsr_r2/$SAMPLE"_recalibration_plot.2.pdf"

    CMD="$SINGULARITY gatk AnalyzeCovariates \
        -before $BEFORE_TABLE \
        -after $AFTER_TABLE \
        -plots $OUT_PDF"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME $LOG $DEPEND"

    echo $CMD >$SCRIPT
    cat $SCRIPT | $JOB_QSUB  

done
