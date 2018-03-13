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
cd $BQSR_DIR

# MERGE ROUND 1 ----------------------------------------------------------------
ROUND=2
mkdir $BQSR_DIR/cohort_vcf_r2

ls genotype_interval_vcfs_r2/*.vcf >cohort_vcf_r2/interval_vcf.list

$SINGULARITY split -d -n l/1000 --additional-suffix .list cohort_vcf_r2/interval_vcf.list cohort_vcf_r2/int.

for INTERVAL in $(seq -w 0 999); do

    MERGE_JOB_NAME="merge_cohort_round_"$INTERVAL"_"$ROUND
    THREADS=1

    IN_LIST="cohort_vcf_r2/int."$INTERVAL".list"
    OUT_VCF="cohort_vcf_r2/$INTERVAL.vcf"

    MERGE_QSUB="$QSUB -pe mpi $THREADS -N $MERGE_JOB_NAME -o logs/$MERGE_JOB_NAME.log"
    MERGE="$SINGULARITY gatk MergeVcfs -I $IN_LIST -O $OUT_VCF -R $REFERENCE"
    echo $MERGE >scripts/$MERGE_JOB_NAME.sh

    cat scripts/$MERGE_JOB_NAME.sh | $MERGE_QSUB

done


#                               <...>
#                       wait till all finished
#                               <...>

#check for completion
PASSED=0
FAILED=0
TOTAL=0
EXPECTED=1000

for INTERVAL in $(seq -w 0 999); do

    MERGE_JOB_NAME="merge_cohort_round_"$INTERVAL"_"$ROUND
    LOG="logs/$MERGE_JOB_NAME.log"
    THREADS=12

    if [[ $(grep "picard.vcf.MergeVcfs done" $LOG) ]]; then
        PASSED=$((PASSED+1))
    else
        FAILED=$((FAILED+1))

        MERGE_QSUB="$QSUB -pe mpi $THREADS -N $MERGE_JOB_NAME -o logs/$MERGE_JOB_NAME.log"

        cat scripts/$MERGE_JOB_NAME.sh | $MERGE_QSUB
    fi

    NUM_SAMPLES=$((NUM_SAMPLES+1))
    TOTAL=$((TOTAL+1))
done

echo -e "PASSED\tFAILED\tTOTAL\tEXPECTED"
echo -e "$PASSED\t$FAILED\t$TOTAL\t$EXPECTED"

#           PASSED  FAILED  TOTAL   EXPECTED
#1stPass    979     21      1000    1000
#2ndPass    1000    0       1000    1000


# MERGE ROUND 2 ----------------------------------------------------------------
ROUND=2

mkdir $BQSR_DIR"/cohort_vcf_r2/round_"$ROUND
cd $BQSR_DIR"/cohort_vcf_r2/round_"$ROUND
mkdir logs scripts

ls $BQSR_DIR/cohort_vcf_r2/*.vcf >interval_vcf.list


$SINGULARITY split -d -n l/100 --additional-suffix .list interval_vcf.list int.

for INTERVAL in $(seq -w 0 99); do

    MERGE_JOB_NAME="merge_cohort_round_"$INTERVAL"_"$ROUND
    THREADS=1

    IN_LIST="int."$INTERVAL".list"
    OUT_VCF="$INTERVAL.vcf"

    MERGE_QSUB="$QSUB -pe mpi $THREADS -N $MERGE_JOB_NAME -o logs/$MERGE_JOB_NAME.log"


    MERGE="$SINGULARITY gatk MergeVcfs -I  $IN_LIST -O $OUT_VCF -R $REFERENCE"
    echo $MERGE >scripts/$MERGE_JOB_NAME.sh

    cat scripts/$MERGE_JOB_NAME.sh | $MERGE_QSUB

done


#                               <...>
#                       wait till all finished
#                               <...>

#check for completion
PASSED=0
FAILED=0
TOTAL=0
EXPECTED=100

for INTERVAL in $(seq -w 0 99); do

    MERGE_JOB_NAME="merge_cohort_round_"$INTERVAL"_2"
    LOG="logs/$MERGE_JOB_NAME.log"
    THREADS=12

    if [[ $(grep "picard.vcf.MergeVcfs done" $LOG) ]]; then
        PASSED=$((PASSED+1))
    else
        FAILED=$((FAILED+1))

        MERGE_QSUB="$QSUB -pe mpi $THREADS -N $MERGE_JOB_NAME -o logs/$MERGE_JOB_NAME.log"

        cat scripts/$MERGE_JOB_NAME.sh | $MERGE_QSUB
    fi

    NUM_SAMPLES=$((NUM_SAMPLES+1))
    TOTAL=$((TOTAL+1))

done

echo -e "PASSED\tFAILED\tTOTAL\tEXPECTED"
echo -e "$PASSED\t$FAILED\t$TOTAL\t$EXPECTED"

#           PASSED  FAILED  TOTAL   EXPECTED
#1stPass    94      6       100     100
#2ndPass    100     0       100     100


# MERGE ROUND 3 ----------------------------------------------------------------
ROUND=3

mkdir $BQSR_DIR"/cohort_vcf_r2/round_"$ROUND
cd $BQSR_DIR"/cohort_vcf_r2/round_"$ROUND

mkdir logs scripts
ls ../round_2/*.vcf >interval_vcf.list

MERGE_JOB_NAME="merge_cohort_round_"$ROUND
THREADS=1

IN_LIST="interval_vcf.list"
OUT_VCF="$BQSR_DIR/cohort_vcf_r2/tmp_cohort.vcf"

MERGE="$SINGULARITY gatk MergeVcfs -I $IN_LIST -O $OUT_VCF -R $REFERENCE"

MERGE_QSUB="$QSUB -pe mpi $THREADS -N $MERGE_JOB_NAME -o logs/$MERGE_JOB_NAME.log"
echo $MERGE >scripts/$MERGE_JOB_NAME.sh

cat scripts/$MERGE_JOB_NAME.sh | $MERGE_QSUB

#                               <...>
#                       wait till all finished
#                               <...>


# SORT COHORT VCF --------------------------------------------------------------
SORT_JOB_NAME=sort_cohort
THREADS=12

cd $BQSR_DIR

IN_GVCF="$BQSR_DIR/cohort_vcf_r2/tmp_cohort.vcf"
OUT_GVCF="$BQSR_DIR/cohort_r1_BQSR.g.vcf"

SORT="$SINGULARITY gatk --java-options "'"-Xmx8G"'" SortVcf -I $IN_GVCF -O $OUT_GVCF"

SORT_QSUB="$QSUB -pe mpi $THREADS -N $SORT_JOB_NAME -o logs/$SORT_JOB_NAME.log"
echo $SORT >scripts/$SORT_JOB_NAME.sh

cat scripts/$SORT_JOB_NAME.sh | $SORT_QSUB
