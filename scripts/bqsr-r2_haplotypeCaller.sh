#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr-r2_haplotypeCaller.sh - use BSQR tables for second round of HC in recalibration proc.

# TODO(neal): update comments


source /master/nplatt/sH_hybridization/scripts/set-env.sh

cd $BQSR_DIR

mkdir hc_vcf_r2

for SAMPLE in $(cat $SAMPLE_LIST); do
    for INTERVAL in $(seq -w 0 49); do    
    
        # HC -------------------------------------------------------------------
        HC_JOB_NAME=$SAMPLE".r2."$INTERVAL.hc
        THREADS=1

        IN_BED="$INTERVALS_DIR/filtered_interval.part$INTERVAL.list"
        IN_BAM=$SAMPLE".bqsr-1.bam"
        OUT_GVCF=$BQSR_DIR"/hc_vcf_r2/"$SAMPLE"_interval_"$INTERVAL".g.vcf"
    
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
        HC_JOB_NAME=$SAMPLE".r2."$INTERVAL.hc
        LOG="logs/$HC_JOB_NAME.log"

        if [[ $(grep "HaplotypeCaller done" $LOG) ]]; then
            PASSED=$((PASSED+1))
        else
            FAILED=$((FAILED+1))

            THREADS=12
            HC_QSUB="$QSUB -pe mpi $THREADS -N $HC_JOB_NAME -o logs/$HC_JOB_NAME.log"
            cat scripts/$HC_JOB_NAME.sh | $HC_QSUB
        fi 

    NUM_SAMPLES=$((NUM_SAMPLES+1))
    TOTAL=$((TOTAL+1))
    done
done 

echo -e "PASSED\tFAILED\tTOTAL\tNUM_SAMPLES\tEXPECTED"
echo -e "$PASSED\t$FAILED\t$TOTAL\t$NUM_SAMPLES\t$EXPECTED"

#        PASSED  FAILED  TOTAL   NUM_SAMPLES     EXPECTED
#1stPass 4799    1       4800    4800            4800
#2ndPass 4800    0       4800    4800            4800


