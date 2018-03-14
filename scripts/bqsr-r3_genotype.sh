#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr-r3_genotype.sh - runs GATK's GenotypeGVCF on each sample/interval 
#   combination

# TODO (nplatt): update comments

#import env variables
source /master/nplatt/sH_hybridization/scripts/set-env.sh

#move to working dir
ROUND=r3
WORK_DIR=$BQSR_DIR/$ROUND"_genotype"

mkdir $WORK_DIR
cd $WORK_DIR

#genotype each interval in the db_r3 dir (qsub)
for INTERVAL in $(ls $BQSR_DIR/$ROUND"_db"); do
    #job specific params    
    JOB_NAME=$INTERVAL".genotype_"$ROUND
    THREADS=1
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND=""
    SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

    #job specific in/out files
    IN_DB=gendb://$BQSR_DIR/$ROUND"_db"/$INTERVAL
    OUT_VCF=$WORK_DIR/$INTERVAL.vcf
    
    #job specific qsub    
    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"
    
    #job specific script
    CMD="$SINGULARITY gatk GenotypeGVCFs \
        -R $REFERENCE \
        -V $IN_DB \
        -new-qual \
        -O $OUT_VCF"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"    

done


#
#                               <...wait...>
#

#sleep while all jobs are running
WAIT_FOR_CLEAR_QUEUE
#
#                               <...wait...>
#


################################################################################
# Check log files to see that all ran to completion 

#reset variables
FAILED="1"
while [ $FAILED -ne 0]; do
    PASSED=0
    FAILED=0
    TOTAL=0
    EXPECTED=$(ls $BQSR_DIR/$ROUND"_db" | wc -l)

    #check that each log file has "GenotypeGVCFs done" indicating run to completion
    for INTERVAL in $(ls $BQSR_DIR/$ROUND"_db"/$INTERVAL); do
        JOB_NAME=$INTERVAL".genotype_"$ROUND
        THREADS=1
        LOG="$LOGS_DIR/$JOB_NAME.log" 
        DEPEND=""
        SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"
        JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

        if [[ $(grep "GenotypeGVCFs done" $LOG) ]]; then
            PASSED=$((PASSED+1))

        else
            #resubmit with 12 threads (only reason doing this is to "hog" memory
            #  from an entire node rather than sharing        
            FAILED=$((FAILED+1))
            THREADS=12
    
            rm $LOG
            cat $SCRIPT | $JOB_QSUB
        fi

        TOTAL=$((TOTAL+1))
    done

    WAIT_FOR_CLEAR_QUEUE

done
