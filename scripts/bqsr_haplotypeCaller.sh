#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr_haplotypeCaller.sh - use BSQR tables for second round of HC in recalibration proc.

# TODO(neal): update comments

cd $BQSR_DIR/$ROUND"_hc"

FAILED="0"
ITERATION="1"

while [ $FAILED -gt 0 ] || [ $ITERATION -le 1 ] ; do

    FAILED=0

    for SAMPLE in $(cat $SAMPLE_LIST); do
        for INTERVAL in $(seq -w 0 49); do    
    
            # HC ---------------------------------------------------------------
            JOB_NAME="snp."$SAMPLE".$ROUND."$INTERVAL.hc
            THREADS=1
            LOG="$LOGS_DIR/$JOB_NAME.log" 
            DEPEND=""
            SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

            IN_BAM=$UNMODIFIED_BAM_DIR/$SAMPLE*".bam"
            IN_BED="$INTERVALS_DIR/filtered_interval.part$INTERVAL.list"
            OUT_GVCF=$BQSR_DIR/$ROUND"_hc"/$SAMPLE"_interval_"$INTERVAL.$ROUND".g.vcf"
    
            JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

            CMD="$SINGULARITY gatk HaplotypeCaller \
                -I $IN_BAM \
                -O $OUT_GVCF \
                -R $REFERENCE \
                -L $IN_BED \
                -ERC GVCF"

        if [ ! -f $LOG ] || ! grep -q "HaplotypeCaller done" $LOG ; then
                FAILED=$((FAILED+1))

                if [ $ITERATION -gt 2 ]; then 
                    THREADS=12;
                fi

                DELETE $LOG $SCRIPT
                SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"

                #only submit a limited number of jobs at a time
                #  ...(dont overload queue)
                LIMIT_RUNNING_JOBS_TO 3000

            fi #closes failed sample loop
        
        done #closes interval loop
    done  #closes sample loop
    
    #wait for queue to clear before checking the number of failed jobs
    WAIT_FOR_CLEAR_QUEUE

    ITERATION=$((ITERATION+1))
done #closes while (failed) loop

