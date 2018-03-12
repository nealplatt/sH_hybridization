#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr-r2_haplotypeCaller.sh - use BSQR tables for second round of HC in recalibration proc.

# TODO(nplatt): update comments
# TODO(nplatt): has not been through test-run yet

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
            #cat scripts/$HC_JOB_NAME.sh | $HC_QSUB
        fi 

    NUM_SAMPLES=$((NUM_SAMPLES+1))
    TOTAL=$((TOTAL+1))
    done
done 

echo -e "PASSED\tFAILED\tTOTAL\tNUM_SAMPLES\tEXPECTED"
echo -e "$PASSED\t$FAILED\t$TOTAL\t$NUM_SAMPLES\t$EXPECTED"

#-------------------------------------------------------------------------------
#NOTES FROM GATK FORUM:
#https://gatkforums.broadinstitute.org/gatk/discussion/4913/recommended-protocol-for-bootstrapping-haplotypecaller-and-baserecalibrator-outputs

#1.1.4 HaplotypeCaller -R $REF -I $BAM_FILE -ERC GVCF -BQSR post_recal_data.1.table -o recal_1.g.vcf
#1.1.5 Filter recal_1.g.vcf --> recal_1.g.filtered.vcf

#Iteration 2
#1.2.1 BaseRecalibrator -R $REF -I $BAM_FILE -knownSites recal_1.g.filtered.vcf -o recal_data.2.table
#1.2.2 BaseRecalibrator -R $REF -I $BAM_FILE -knownSites recal_1.g.filtered.vcf -BQSR recal_data.2.table -o post_recal_data.2.table
#1.2.3 AnalyzeCovariates -R $REF -before recal_data.2.table -after post_recal_data.2.table -plots recalibration_plots.2.pdf
#1.2.4 HaplotypeCaller -R $REF -I $BAM_FILE -ERC GVCF -BQSR post_recal_data.2.table -o recal_2.g.vcf
#1.2.5 Filter recal_2.g.vcf --> recal_2.g.filtered.vcf

#--> HC -BQSR on the fly recalibration of "raw" bam files
#--> Each iteration should produce a new vcf file that is used for the recalibration of "raw" reads
#--> Differences between the AnalyzeCovariates before and after should diminush after each iteration (difference 1.1.3 > difference 1.2.3)
