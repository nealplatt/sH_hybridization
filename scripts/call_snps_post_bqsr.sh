#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# haplotypeCaller.sh - use recalibrated reads to call SNPS.

# TODO(neal): update comments

source /master/nplatt/sH_hybridization/scripts/set-env.sh

RECALIBRATED_BAM_DIR=$BQSR_DIR/r1_bqsr_bams

mkdir $RESULTS_DIR/haplotype_caller
cd $RESULTS_DIR/haplotype_caller


# Haplotype Caller -------------------------------------------------------------
FAILED="0"
ITERATION="1"
while [ $FAILED -gt 0 ] || [ $ITERATION -le 1 ] ; do

    FAILED=0

    for SAMPLE in $(cat $SAMPLE_LIST); do
        for INTERVAL in $(seq -w 0 49); do    
    
            # HC ---------------------------------------------------------------
            JOB_NAME="snp."$SAMPLE".final."$INTERVAL.hc
            THREADS=1
            LOG="$LOGS_DIR/$JOB_NAME.log" 
            DEPEND=""
            SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

            IN_BAM=$RECALIBRATED_BAM_DIR/$SAMPLE*".bam"
            IN_BED="$INTERVALS_DIR/filtered_interval.part$INTERVAL.list"
            OUT_GVCF=$SAMPLE"_interval_"$INTERVAL."final.g.vcf"
    
            JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

            CMD="$SINGULARITY gatk HaplotypeCaller \
                -I $IN_BAM \
                -O $OUT_GVCF \
                -R $REFERENCE \
                -L $IN_BED \
                -ERC GVCF"

            if [ ! -f $LOG ] || ! grep -q "HaplotypeCaller done" $LOG; then
                FAILED=$((FAILED+1))
                if [ $ITERATION -gt 2 ]; then 
                    THREADS=12;
                fi
                
                DELETE $LOG $SCRIPT
                SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"
        
        done #closes interval loop
    done  #closes sample loop
    
    #wait for queue to clear before checking the number of failed jobs
    WAIT_FOR_CLEAR_QUEUE

    ITERATION=$((ITERATION+1))
done #closes while (failed) loop



# Make indivudual vcfs ---------------------------------------------------------
FAILED="1"
while [ $FAILED -ne 0 ]; do
    FAILED=0

    for SAMPLE in $(cat $SAMPLE_LIST); do

            # MERGE_INDIV_GVCF -------------------------------------------------
            ls  $SAMPLE"_interval_"*".g.vcf" >$SAMPLE.list      

            JOB_NAME="snp."$SAMPLE.final.merge
            THREADS=12
            LOG="$LOGS_DIR/$JOB_NAME.log" 
            DEPEND=""
            SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

            IN_LIST=$SAMPLE.list
            OUT_GVCF="tmp_"$SAMPLE".merged.g.vcf"
    
            JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

            CMD="$SINGULARITY gatk MergeVcfs -I $IN_LIST -O $OUT_GVCF"

            if [ ! -f $LOG ] || ! grep -q "picard.vcf.MergeVcfs done" $LOG; then
                FAILED=$((FAILED+1))
                DELETE $LOG $SCRIPT
                SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"
            fi

            # SORT_INDIV_GVCF --------------------------------------------------
            JOB_NAME="snp.$SAMPLE.final.sort"
            THREADS=12
            LOG="$LOGS_DIR/$JOB_NAME.log" 
            DEPEND="-hold_jid snp.$SAMPLE.final.merge"
            SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

            IN_GVCF=$OUT_GVCF
            OUT_GVCF=$SAMPLE"_final.g.vcf"
    
            JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

            CMD="$SINGULARITY gatk SortVcf -I $IN_GVCF -O $OUT_GVCF"
            
            if [ ! -f $LOG ] || ! grep -q "picard.vcf.SortVcf done" $LOG; then
                FAILED=$((FAILED+1))
                DELETE $LOG $SCRIPT
                SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"
            fi
    
    done
    #sleep while all jobs are running   
    WAIT_FOR_CLEAR_QUEUE

done



#clean up real quick
rm tmp_*.merged.g.vcf* *.list

mkdir hc_vcfs
mv *_interval_??.final.g.vcf* hc_vcfs
tar -czf hc_vcfs.tgz hc_vcfs
rm -r hc_vcfs



# GDBIMPORT --------------------------------------------------------------------
ls $(pwd)/*"_final.g.vcf" >final.gvcf.list
mkdir db

# loop for submission
FAILED="1"
while [ $FAILED -ne 0 ]; do
    FAILED=0
    
    for INTERVAL in $(cat $INTERVALS_DIR/all_filtered_intervals.list); do

        SAFE_INTERVAL_NAME=$(echo $INTERVAL | sed 's/:/-/')

        JOB_NAME="snp."$SAFE_INTERVAL_NAME"_final"
        THREADS=12
        LOG="$LOGS_DIR/$JOB_NAME.log" 
        DEPEND=""
        SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

        IN_LIST=final.gvcf.list 
        OUT_DB=$(pwd)/"db/$SAFE_INTERVAL_NAME"
    
        JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

        CMD="$SINGULARITY gatk --java-options "'"-Xmx4g -Xms4g"'" \
                GenomicsDBImport \
                -V $IN_LIST \
                --genomicsdb-workspace-path $OUT_DB \
                -L $INTERVAL \
                --reader-threads $THREADS \
                --batch-size 24"

        if [ ! -f $LOG ] || ! grep -q "genomicsdb.GenomicsDBImport done" $LOG; then

            LIMIT_RUNNING_JOBS_TO 300

            FAILED=$((FAILED+1))
            DELETE $LOG $SCRIPT
            SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"
        fi
    done
    
    #sleep while all jobs are running   
    WAIT_FOR_CLEAR_QUEUE

done

# GENOTYPE ---------------------------------------------------------------------
FAILED="1"
while [ $FAILED -ne 0 ]; do
    FAILED=0

    for INTERVAL in $(ls db); do

    #job specific params    
    JOB_NAME="snp."$INTERVAL".genotype_"$ROUND
    THREADS=1
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND=""
    SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

    #job specific in/out files
    IN_DB=gendb://$(pwd)/db/$INTERVAL
    OUT_VCF=$INTERVAL.vcf
    
    #job specific qsub    
    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"
    
    #job specific script
    CMD="$SINGULARITY gatk GenotypeGVCFs \
        -R $REFERENCE \
        -V $IN_DB \
        -new-qual \
        -O $OUT_VCF"

        if [ ! -f $LOG ] || ! grep -q "GenotypeGVCFs done" $LOG; then

            FAILED=$((FAILED+1))
            DELETE $LOG $SCRIPT
            SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"
        fi  
    done
    
    #sleep while all jobs are running   
    WAIT_FOR_CLEAR_QUEUE

done


#after this has successfully run will need to 
# 1) build a cohort vcf
# 2) filter/refine genotypes

#When HC is finished 
# 1) filter SNPS from individual vcfs
