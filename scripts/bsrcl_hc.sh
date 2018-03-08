#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bsrcl_hc.sh - iterates haplotype caller across all samples and individuals

# TODO(nplatt): Update comments

#load variables
source master/nplatt/sH_hybridizationscripts/set_env.sh

#Filter reads wit trimmomatic
mkdir $BSRCL_DIR $BSRCL_DIR/scripts $BSRCL_DIR/logs

cd $BSRCL_DIR

for SAMPLE in $(cat $SAMPLE_LIST); do
    for INTERVAL in $(seq -w 0 49); do    
    
        # HC -------------------------------------------------------------------
        HC_JOB_NAME=$SAMPLE"."$INTERVAL.hc
        THREADS=1

        IN_BED="$INTERVALS_DIR/filtered_interval.part$INTERVAL.list"
        IN_BAM=$MAP_DIR/$SAMPLE"_processed.bam"
        OUT_GVCF=$BSRCL_DIR/$SAMPLE"_interval_"$INTERVAL".g.vcf"
    
        HC="$SINGULARITY gatk HaplotypeCaller \
            -I $IN_BAM \
            -O $OUT_GVCF \
            -R $REFERENCE \
            -L $IN_BED \
            -ERC GVCF"

        HC_QSUB="$QSUB -pe mpi $THREADS -N $HC_JOB_NAME -o logs/$HC_JOB_NAME.log"
        echo $HC >scripts/$HC_JOB_NAME.sh
    
        cat scripts/$HC_JOB_NAME.sh | $HC_QSUB

    done
done


# HC-CHECKER------------------------------------------------------------
PASSED=0
FAILED=0
TOTAL=0
NUM_SAMPLES=0
EXPECTED=$(expr 96 \* 50)

for SAMPLE in $(cat $SAMPLE_LIST); do
    
    for INTERVAL in $(seq -w 0 49); do
        HC_JOB_NAME=$SAMPLE"."$INTERVAL.hc
        LOG="logs/$HC_JOB_NAME.log"

        if [[ $(grep "HaplotypeCaller done" $LOG) ]]; then
            PASSED=$((PASSED+1))
        else
            FAILED=$((FAILED+1))
            HC_JOB_NAME=$SAMPLE"."$INTERVAL.hc
            THREADS=4
            HC_QSUB="$QSUB -pe mpi $THREADS -N $HC_JOB_NAME -o logs/$HC_JOB_NAME.log"
            cat scripts/$HC_JOB_NAME.sh | $HC_QSUB
        fi 

    NUM_SAMPLES=$((NUM_SAMPLES+1))
    TOTAL=$((TOTAL+1))
    done
done 

echo -e "PASSED\tFAILED\tTOTAL\tNUM_SAMPLES\tEXPECTED"
echo -e "$PASSED\t$FAILED\t$TOTAL\t$NUM_SAMPLES\t$EXPECTED"

#re-run untill all samples completed
