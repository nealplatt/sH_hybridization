#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr-r2_genotype.sh - runs GATK's GenotypeGVCF on each sample/interval 
#   combination

# TODO(nplatt): update comments

#import env variables
source /master/nplatt/sH_hybridization/scripts/set-env.sh

#move to working dir
cd $BQSR_DIR

#all output will go to this dir
mkdir genotype_interval_vcfs_r2

#genotype each interval in the db_r2 dir (qsub)
for INTERVAL in $(ls db_r2); do
    #job specific params    
    GENOTYPE_JOB_NAME=$INTERVAL".genotype_r2"
    THREADS=1
    LOG="-o logs/$GENOTYPE_JOB_NAME.log"
    DEPEND=""

    #job specific in/out files
    IN_DB="gendb://db_r2/$INTERVAL"
    OUT_VCF="genotype_interval_vcfs_r2/$INTERVAL.vcf"
    
    #job specific script
    GENOTYPE="$SINGULARITY gatk GenotypeGVCFs \
        -R $REFERENCE \
        -V $IN_DB \
        -new-qual \
        -O $OUT_VCF"
    echo $GENOTYPE >scripts/$GENOTYPE_JOB_NAME.sh
    
    #job specific qsub
    JOB_QSUB="$QSUB -pe mpi $THREADS -N $GENOTYPE_JOB_NAME $LOG $DEPEND"

    #submit job to scheduler
    cat scripts/$GENOTYPE_JOB_NAME.sh | $JOB_QSUB
    
done

################################################################################
# Check log files to see that all ran to completion 

#reset variables
PASSED=0
FAILED=0
TOTAL=0
EXPECTED=$(ls -d db_r2 | wc -l)

#check that each log file has "GenotypeGVCFs done" indicating run to completion
for INTERVAL in $(ls db_r2); do
    GENOTYPE_JOB_NAME=$INTERVAL".genotype_r2"
    THREADS=1
    LOG="-o logs/$GENOTYPE_JOB_NAME.log"
    DEPEND=""


    if [[ $(grep "GenotypeGVCFs done" $LOG) ]]; then
        PASSED=$((PASSED+1))
    else
        #resubmit with 12 threads (only reason doing this is to "hog" memory
        #  from an entire node rather than sharing        
        FAILED=$((FAILED+1))
        THREADS=12
    
        JOB_QSUB="$QSUB -pe mpi $THREADS -N $GENOTYPE_JOB_NAME $LOG $DEPEND"
        cat scripts/$GENOTYPE_JOB_NAME.sh | $GENOTYPE_QSUB
    fi

    TOTAL=$((TOTAL+1))
done

echo -e "PASSED\tFAILED\tTOTAL\tEXPECTED"
echo -e "$PASSED\t$FAILED\t$TOTAL\t$EXPECTED"

#           PASSED  FAILED  TOTAL   EXPECTED
#1stPass

