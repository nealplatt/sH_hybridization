#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr-r3_prepForGenotype.sh - takes VCF files from haplotype caller and preps 
#   for genotyping by creating 1 vcf per individual and importing into a DB

# TODO(nplatt): first pass

source /master/nplatt/sH_hybridization/scripts/set-env.sh

ROUND=r3

WORK_DIR=$BQSR_DIR/$ROUND"_indvi_vcf"

mkdir $WORK_DIR
cd $WORK_DIR

for SAMPLE in $(cat $SAMPLE_LIST); do

    # MERGE_INDIV_GVCF ---------------------------------------------------------
    ls  $BQSR_DIR/r3_hc/$SAMPLE"_interval_"*".g.vcf" >$SAMPLE.list      

    JOB_NAME=$SAMPLE.$ROUND.merge
    THREADS=12
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND=""
    SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

    IN_LIST=$SAMPLE.list
    OUT_GVCF="tmp_"$SAMPLE".merged.g.vcf"
    
    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    CMD="$SINGULARITY gatk MergeVcfs -I $IN_LIST -O $OUT_GVCF"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"


    # SORT_INDIV_GVCF ----------------------------------------------------------
    JOB_NAME=$SAMPLE.$ROUND.sort
    THREADS=12
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid $SAMPLE.$ROUND.merge"
    SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

    IN_GVCF=$OUT_GVCF
    OUT_GVCF=$SAMPLE"_bqsr-"$ROUND".g.vcf"
    
    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    CMD="$SINGULARITY gatk SortVcf -I $IN_GVCF -O $OUT_GVCF"

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


#check that sort and merge steps finished to completion
FAILED="1"
while [ $FAILED -ne 0]; do
    PASSED=0
    FAILED=0
    TOTAL=0
    EXPECTED=96

    for SAMPLE in $(cat $SAMPLE_LIST); do
        THREADS=12
                
        MERGE_JOB_NAME=$SAMPLE.$ROUND.merge
        SORT_JOB_NAME=$SAMPLE.$ROUND.sort
      
        MERGE_LOG="$LOGS_DIR/$MERGE_JOB_NAME.log"  
        SORT_LOG="$LOGS_DIR/$SORT_JOB_NAME.log"

        MERGE_SCRIPT="$SCRIPTS_DIR/$MERGE_JOB_NAME.sh"  
        SORT_SCRIPT="$SCRIPTS_DIR/$SORT_JOB_NAME.sh"

        MERGE_QSUB="$QSUB -pe mpi $THREADS -N $MERGE_JOB_NAME -o $MERGE_LOG"
        SORT_QSUB="$QSUB -pe mpi $THREADS -N $SORT_JOB_NAME -o $SORT_LOG -hold_jid $MERGE_JOB_NAME"


        if [[ $(grep "picard.vcf.MergeVcfs done" $MERGE_LOG) && $(grep "picard.vcf.SortVcf done" $SORT_LOG) ]]; then
            PASSED=$((PASSED+1))
        else
            FAILED=$((FAILED+1))

            rm $MERGE_LOG $SORT_LOG
        
            cat $MERGE_SCRIPT | $MERGE_QSUB
            cat $SORT_SCRIPT | $SORT_QSUB
        fi 

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


# GDBIMPORT ----------------------------------------------------------------
DB_DIR=$BQSR_DIR/$ROUND"_db"

mkdir $DB_DIR

for SAMPLE in $(cat $SAMPLE_LIST); do
    echo $BQSR_DIR/$ROUND"_indvi_vcf"/$SAMPLE"_bqsr-"$ROUND".g.vcf" >>$WORK_DIR/samples_$ROUND.list
done

# loop for submission
for INTERVAL in $(cat $INTERVALS_DIR/all_filtered_intervals.list); do

    SAFE_INTERVAL_NAME=$(echo $INTERVAL | sed 's/:/-/')

    JOB_NAME=$SAFE_INTERVAL_NAME"_"$ROUND
    THREADS=12
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND=""
    SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

    IN_LIST=$WORK_DIR/samples_$ROUND.list  
    OUT_DB="$DB_DIR/$SAFE_INTERVAL_NAME"
    
    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    CMD="$SINGULARITY gatk --java-options "'"-Xmx4g -Xms4g"'" \
        GenomicsDBImport \
        -V $WORK_DIR/samples_$ROUND.list \
        --genomicsdb-workspace-path $OUT_DB \
        -L $INTERVAL \
        --reader-threads $THREADS \
        --batch-size 24"

    DELETE $LOG $SCRIPT

    #only submit a limited number of jobs at a time...(dont overload queue)
    LIMIT_RUNNING_JOBS_TO 300
        
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


#CHECK FOR COMPLETION OF IMPORT
FAILED="1"
while [ $FAILED -ne 0]; do
    PASSED=0
    FAILED=0
    TOTAL=0
    EXPECTED=$(wc -l $INTERVALS_DIR/all_filtered_intervals.list)


    for INTERVAL in $(cat $INTERVALS_DIR/all_filtered_intervals.list); do
        SAFE_INTERVAL_NAME=$(echo $INTERVAL | sed 's/:/-/')
        JOB_NAME=$SAFE_INTERVAL_NAME"_"$ROUND
        THREADS=12
        LOG="$LOGS_DIR/$JOB_NAME.log" 
        DEPEND=""
        SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"
        JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

        if [[ $(grep "genomicsdb.GenomicsDBImport done" $LOG) ]]; then
            PASSED=$((PASSED+1))
        else
            FAILED=$((FAILED+1))
       
            rm -r $DB_DIR/$SAFE_INTERVAL_NAME 
            rm $LOG
            cat $SCRIPT | $JOB_QSUB

        fi 

        TOTAL=$((TOTAL+1))
    done
done

