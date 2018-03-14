#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr-r2_buildCohortVcf.sh - after genotyping, build and filter cohort vcf for 
#       recalibration

# TODO(nplatt): update comments

# done manually/interactivley on the head/scheduler node --- sue me.
source /master/nplatt/sH_hybridization/scripts/set-env.sh

ROUND=r3
WORK_DIR=$BQSR_DIR/$ROUND"_bqsr"

mkdir $WORK_DIR
cd $WORK_DIR


# MERGE ROUND 1 ----------------------------------------------------------------
ls $BQSR_DIR/$ROUND"_genotype"/*.vcf >interval_vcf.list

$SINGULARITY split -d -n l/1000 --additional-suffix .list interval_vcf.list int.

for INTERVAL in $(seq -w 0 999); do

    JOB_NAME="merge_cohort_"$INTERVAL"_"$ROUND
    THREADS=1
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND=""
    SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

    IN_LIST="$WORK_DIR/int."$INTERVAL".list"
    OUT_VCF="$WORK_DIR/$INTERVAL.vcf"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    CMD="$SINGULARITY gatk MergeVcfs -I $IN_LIST -O $OUT_VCF -R $REFERENCE"

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
FAILED="1"
while [ $FAILED -ne 0 ]; do
    PASSED=0
    FAILED=0
    TOTAL=0
    EXPECTED=1000

    for INTERVAL in $(seq -w 0 999); do
        JOB_NAME="merge_cohort_"$INTERVAL"_"$ROUND
        THREADS=1
        LOG="$LOGS_DIR/$JOB_NAME.log" 
        DEPEND=""
        SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"
        JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

        if [[ $(grep "picard.vcf.MergeVcfs done" $LOG) ]]; then
            PASSED=$((PASSED+1))
        else
            #resubmit with 12 threads (only reason doing this is to "hog" memory
            #  from an entire node rather than sharing        
            FAILED=$((FAILED+1))
            THREADS=12
    
            rm $LOG
            cat $SCRIPT | $JOB_QSUB
        fi

        NUM_SAMPLES=$((NUM_SAMPLES+1))
        TOTAL=$((TOTAL+1))
    done
done

#
#                               <...wait...>
#
#sleep while all jobs are running
WAIT_FOR_CLEAR_QUEUE
#
#                               <...wait...>
#

# MERGE ROUND 2 ----------------------------------------------------------------
rm *.list
ls $WORK_DIR/*.vcf >interval_vcf.list

$SINGULARITY split -d -n l/100 --additional-suffix .list interval_vcf.list int2.

for INTERVAL in $(seq -w 0 99); do

    JOB_NAME="merge2_cohort_"$INTERVAL"_"$ROUND
    THREADS=1
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND=""
    SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

    IN_LIST="$WORK_DIR/int2."$INTERVAL".list"
    OUT_VCF="$WORK_DIR/$INTERVAL.2.vcf"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    CMD="$SINGULARITY gatk MergeVcfs -I $IN_LIST -O $OUT_VCF -R $REFERENCE"

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


# MERGE ROUND 3 ----------------------------------------------------------------
rm *.list
ls $WORK_DIR/*.2.vcf >interval3_vcf.list

JOB_NAME="merge3_cohort_"$INTERVAL"_"$ROUND
THREADS=1
LOG="$LOGS_DIR/$JOB_NAME.log" 
DEPEND=""
SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

IN_LIST="$WORK_DIR/interval3_vcf.list"
OUT_VCF="$WORK_DIR/tmp_cohort.vcf"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

CMD="$SINGULARITY gatk MergeVcfs -I $IN_LIST -O $OUT_VCF -R $REFERENCE"

DELETE $LOG $SCRIPT
SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"

# SORT COHORT VCF --------------------------------------------------------------
SORT_JOB_NAME=sort_cohort
THREADS=12
LOG="$LOGS_DIR/$JOB_NAME.log"
DEPEND="-hold_jid merge3_cohort_"$INTERVAL"_"$ROUND
SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

IN_GVCF="$WORK_DIR/tmp_cohort.vcf"
OUT_GVCF="$WORK_DIR/cohort_$ROUND.g.vcf"

CMD="$SINGULARITY gatk --java-options "'"-Xmx8G"'" SortVcf -I $IN_GVCF -O $OUT_GVCF"

SORT_QSUB="$QSUB -pe mpi $THREADS -N $SORT_JOB_NAME -o logs/$SORT_JOB_NAME.log"

DELETE $LOG $SCRIPT
SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"


