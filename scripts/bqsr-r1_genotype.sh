#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bsrcl_genotype - runs GATK's GenotypeGVCF on each sample/interval combination

# TODO(nplatt): add code to check that all combos run successfully
# TODO(nplatt): update comments


source /master/nplatt/sH_hybridization/scripts/set-env.sh

cd $BQSR_DIR

NUM_JOBS_IN_QUEUE=0
MAX_JOBS_ALLOWED=3000

mkdir genotype_interval_vcfs

for INTERVAL in $(ls db); do
    GENOTYPE_JOB_NAME=$INTERVAL".genotype"
    THREADS=1

    IN_DB="gendb://db/$INTERVAL"
    OUT_VCF="genotype_interval_vcfs/$INTERVAL.vcf"
    
    GENOTYPE="$SINGULARITY gatk GenotypeGVCFs -R $REFERENCE -V $IN_DB -new-qual -O $OUT_VCF"

    GENOTYPE_QSUB="$QSUB -pe mpi $THREADS -N $GENOTYPE_JOB_NAME -o logs/$GENOTYPE_JOB_NAME.log"
    echo $GENOTYPE >scripts/$GENOTYPE_JOB_NAME.sh


    NUM_JOBS_IN_QUEUE=$(qstat | grep nplatt | wc -l)

    while [ $NUM_JOBS_IN_QUEUE -gt $MAX_JOBS_ALLOWED ]; do
        sleep 1s
        echo -n "."
        NUM_JOBS_IN_QUEUE=$(qstat | grep nplatt | wc -l)
    done

    cat scripts/$GENOTYPE_JOB_NAME.sh | $GENOTYPE_QSUB
    
done


#################################

#CHECK FOR COMPLETION

PASSED=0
FAILED=0
TOTAL=0
NUM_SAMPLES=0
EXPECTED=$(expr 96 \* 50)

for INTERVAL in $(ls db); do
    GENOTYPE_JOB_NAME=$INTERVAL".genotype"
    LOG="logs/$GENOTYPE_JOB_NAME.log"


    if [[ $(grep "GenotypeGVCFs done" $LOG) ]]; then
        PASSED=$((PASSED+1))
    else
        FAILED=$((FAILED+1))
        THREADS=4
    
        GENOTYPE_QSUB="$QSUB -pe mpi $THREADS -N $GENOTYPE_JOB_NAME -o logs/$GENOTYPE_JOB_NAME.log"
        cat scripts/$GENOTYPE_JOB_NAME.sh | $GENOTYPE_QSUB
    fi

    NUM_SAMPLES=$((NUM_SAMPLES+1))
    TOTAL=$((TOTAL+1))
done

echo -e "PASSED\tFAILED\tTOTAL\tNUM_SAMPLES\tEXPECTED"
echo -e "$PASSED\t$FAILED\t$TOTAL\t$NUM_SAMPLES\t$EXPECTED"
