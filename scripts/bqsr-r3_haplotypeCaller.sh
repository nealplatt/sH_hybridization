#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr-r3_haplotypeCaller.sh - use BSQR tables for second round of HC in recalibration proc.

# TODO(neal): update comments


source /master/nplatt/sH_hybridization/scripts/set-env.sh

ROUND=r3

WORK_DIR=$BQSR_DIR/$ROUND"_hc"

mkdir $WORK_DIR
cd $WORK_DIR

for SAMPLE in $(cat $SAMPLE_LIST); do
    for INTERVAL in $(seq -w 0 49); do    
    
        # HC -------------------------------------------------------------------
        JOB_NAME=$SAMPLE".$ROUND."$INTERVAL.hc
        THREADS=1
        LOG="$LOGS_DIR/$JOB_NAME.log" 
        DEPEND=""
        SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

        IN_BED="$INTERVALS_DIR/filtered_interval.part$INTERVAL.list"
        IN_BAM=$BQSR_DIR/r2_bqsr_bams/$SAMPLE".bqsr-2.bam"
        OUT_GVCF=$WORK_DIR/$SAMPLE"_interval_"$INTERVAL.$ROUND".g.vcf"
    
        JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

        CMD="$SINGULARITY gatk HaplotypeCaller \
            -I $IN_BAM \
            -O $OUT_GVCF \
            -R $REFERENCE \
            -L $IN_BED \
            -ERC GVCF"

        DELETE $LOG $SCRIPT
        SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"

    done
done

#
#                               <...wait...>
#

#sleep while all jobs are running
NUM_JOBS_IN_QUEUE=$(qstat | grep Sh. | wc -l)

   while [ $NUM_JOBS_IN_QUEUE -gt 0 ]; do
        sleep 60s
        echo -n "."
        NUM_JOBS_IN_QUEUE=$(qstat | grep nplatt | wc -l)
    done

#
#                               <...wait...>
#

# HC-CHECKER------------------------------------------------------------
FAILED="1"
while [ $FAILED -ne 0]; do
    PASSED=0
    FAILED=0
    TOTAL=0
    NUM_SAMPLES=0
    EXPECTED=$(expr 96 \* 50)

    for SAMPLE in $(cat $SAMPLE_LIST); do
    
        for INTERVAL in $(seq -w 0 49); do
            JOB_NAME=$SAMPLE".$ROUND."$INTERVAL.hc
            LOG="$LOGS_DIR/$JOB_NAME.log" 

            if [[ $(grep "HaplotypeCaller done" $LOG) ]]; then
                PASSED=$((PASSED+1))
            else
                FAILED=$((FAILED+1))

                THREADS=12
                $JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG"
            
                DELETE $LOG $SCRIPT
                SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"
            fi 

            NUM_SAMPLES=$((NUM_SAMPLES+1))
            TOTAL=$((TOTAL+1))
        done
    done 
done

echo -e "PASSED\tFAILED\tTOTAL\tNUM_SAMPLES\tEXPECTED"
echo -e "$PASSED\t$FAILED\t$TOTAL\t$NUM_SAMPLES\t$EXPECTED"

#        PASSED  FAILED  TOTAL   NUM_SAMPLES     EXPECTED
#1stPass 
#2ndPass 


