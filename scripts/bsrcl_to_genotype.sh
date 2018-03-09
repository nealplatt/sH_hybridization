#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bsrcl_to_genotype.sh - takes VCF files from haplotype caller and preps for 
#  genotyping by creating 1 vcf per individual and importing into a DB

# TODO(nplatt): write code to check log files for successful GDBimport


source master/nplatt/sH_hybridizationscripts/set_env.sh

cd $BSRCL_DIR


MAX_JOBS_ALLOWED=300
NUM_JOBS_IN_QUEUE=$(qstat | grep nplatt | wc -l)

for SAMPLE in $(cat $SAMPLE_LIST); do

    # MERGE_INDIV_GVCF ---------------------------------------------------------
    ls $SAMPLE"_interval_"*".g.vcf" >$SAMPLE.list      

    MERGE_JOB_NAME=$SAMPLE".merge_indiv"
    THREADS=12

    IN_LIST=$SAMPLE.list
    OUT_GVCF="tmp_$SAMPLE.merged.g.vcf"
    
    MERGE="$SINGULARITY gatk MergeVcfs -I $IN_LIST -O $OUT_GVCF"

    MERGE_QSUB="$QSUB -pe mpi $THREADS -N $MERGE_JOB_NAME -o logs/$MERGE_JOB_NAME.log"
    echo $MERGE >scripts/$MERGE_JOB_NAME.sh

    cat scripts/$MERGE_JOB_NAME.sh | $MERGE_QSUB


    # SORT_INDIV_GVCF ----------------------------------------------------------
    SORT_JOB_NAME=$SAMPLE".sort_indiv"
    THREADS=12

    IN_GVCF=$OUT_GVCF
    OUT_GVCF="$SAMPLE.g.vcf"
    
    SORT="$SINGULARITY gatk SortVcf -I $IN_GVCF -O $OUT_GVCF"

    SORT_QSUB="$QSUB -pe mpi $THREADS -N $SORT_JOB_NAME -o logs/$SORT_JOB_NAME.log -hold_jid $MERGE_JOB_NAME"
    echo $SORT >scripts/$SORT_JOB_NAME.sh

    cat scripts/$SORT_JOB_NAME.sh | $SORT_QSUB

done


#CHECK FOR COMPLETION OF SORT AND MERGE STEPS
PASSED=0
FAILED=0
TOTAL=0
EXPECTED=96


for SAMPLE in $(cat $SAMPLE_LIST); do
    MERGE_JOB_NAME=$SAMPLE".merge_indiv"
    SORT_JOB_NAME=$SAMPLE".sort_indiv"
      
    MERGE_LOG="logs/$MERGE_JOB_NAME.log"  
    SORT_LOG="logs/$SORT_JOB_NAME.log"


    if [[ $(grep "picard.vcf.MergeVcfs done" $MERGE_LOG) && $(grep "picard.vcf.SortVcf done" $SORT_LOG) ]]; then
        PASSED=$((PASSED+1))
    else
        FAILED=$((FAILED+1))

        THREADS=12

        MERGE_JOB_NAME=$SAMPLE".merge_indiv"
        SORT_JOB_NAME=$SAMPLE".sort_indiv"

  
        MERGE_QSUB="$QSUB -pe mpi $THREADS -N $MERGE_JOB_NAME -o logs/$MERGE_JOB_NAME.log"
        SORT_QSUB="$QSUB -pe mpi $THREADS -N $SORT_JOB_NAME -o logs/$SORT_JOB_NAME.log -hold_jid $MERGE_JOB_NAME"
        
        cat scripts/$MERGE_JOB_NAME.sh | $MERGE_QSUB
        cat scripts/$SORT_JOB_NAME.sh | $SORT_QSUB
    fi 

    TOTAL=$((TOTAL+1))
done

echo -e "PASSED\tFAILED\tTOTAL\tEXPECTED"
echo -e "$PASSED\t$FAILED\t$TOTAL\t$EXPECTED"


# GDBIMPORT ----------------------------------------------------------------
mkdir db

for INTERVAL in $(cat $INTERVALS_DIR/all_filtered_intervals.list); do
    GDBIMPORT_JOB_NAME=$(echo $SAMPLE.$INTERVAL | sed 's/:/-/')
    THREADS=12

    IN_GVCF="$SAMPLE.g.vcf"
    OUT_DB="db/$INTERVAL"
    

    
    GDBIMPORT="$SINGULARITY gatk --java-options "'"-Xmx4g -Xms4g"'" \
        GenomicsDBImport \
        -V $SAMPLE_LIST \
        --genomicsdb-workspace-path $OUT_DB \
        -L $INTERVAL \
        --reader-threads $THREADS \
        --batch-size 24"

    GDBIMPORT_QSUB="$QSUB -pe mpi $THREADS -N $GDBIMPORT_JOB_NAME -o logs/$GDBIMPORT_JOB_NAME -hold_jid $SORT_JOB_NAME"
    echo $GDBIMPORT >scripts/$GDBIMPORT_JOB_NAME.sh

    #only submit a limited number of jobs at a time...(dont overload queue)
    NUM_JOBS_IN_QUEUE=$(qstat | grep nplatt | wc -l)

    while [ $NUM_JOBS_IN_QUEUE -gt $MAX_JOBS_ALLOWED ]; do
        sleep 1s
        echo -n "."
        NUM_JOBS_IN_QUEUE=$(qstat | grep nplatt | wc -l)
    done
        
    cat scripts/$GDBIMPORT_JOB_NAME.sh | $GDBIMPORT_QSUB
done

#CHECK FOR COMPLETION


