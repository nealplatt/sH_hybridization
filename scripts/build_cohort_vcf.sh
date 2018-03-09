#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# build_cohort_vcf - after genotyping, build and filter cohort vcf for 
#       recalibration

# TODO(nplatt): fix entire script
# TODO(nplatt): update comments

# done manually/interactivley on the head node --- sue me.
source master/nplatt/sH_hybridizationscripts/set_env.sh
cd $BSRCL_DIR

# MERGE ROUND 1 ----------------------------------------------------------------
ROUND=1
mkdir $BSRCL_DIR/round_$ROUND
cd $BSRCL_DIR/round_$ROUND
ls interval_vcf/*.vcf >interval_vcf.list

$SINGULARITY split -d -n l/1000 --additional-suffix .intList interval_vcf.list



for INTERVAL in $(seq -w 1 1000); do

    MERGE_JOB_NAME=$SAMPLE".merge_cohort_round_"$ROUND
    THREADS=1

    IN_LIST="$INTERVAL.intList"
    OUT_VCF="round_$ROUND/$INTERVAL.vcf"
    
    MERGE="gatk --java-options "-Xmx4G" \
    MergeVcfs \
    -I  $IN_LIST \
    -O $OUT_VCF \
    -R $REFERENCE"

    MERGE="gatk --java-options "'"-Xmx4G"'" MergeVcfs -I  $IN_LIST -O $OUT_VCF -R $REFERENCE"
    echo $MERGE >scripts/$MERGE_JOB_NAME.sh

    cat scripts/$MERGE_JOB_NAME.sh | $MERGE_QSUB
 
done


#<<<<<<<<<<<<<<<<<<< wait till all finished >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# MERGE ROUND 2 ----------------------------------------------------------------
ROUND=2
mkdir $BSRCL_DIR/round_$ROUND
cd $BSRCL_DIR/round_$ROUND

ls $BSRCL_DIR/round_1/*.vcf >interval_vcf.list

$SINGULARITY split -d -n l/250 --additional-suffix .intList interval_vcf.list

for INTERVAL in $(seq -w 1 250); do

    MERGE_JOB_NAME=$SAMPLE".merge_cohort_round_"$ROUND
    THREADS=1

    IN_LIST="$INTERVAL.intList"
    OUT_VCF="round_$ROUND/$INTERVAL.vcf"
    
    MERGE="gatk --java-options "'"-Xmx4G"'" MergeVcfs -I  $IN_LIST -O $OUT_VCF -R $REFERENCE"

    MERGE_QSUB="$QSUB -pe mpi $THREADS -N $MERGE_JOB_NAME -o logs/$MERGE_JOB_NAME.log"
    echo $MERGE >scripts/$MERGE_JOB_NAME.sh

    cat scripts/$MERGE_JOB_NAME.sh | $MERGE_QSUB
 
done
#<<<<<<<<<<<<<<<<<<< wait till all finished >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# MERGE ROUND 3 ----------------------------------------------------------------
ROUND=3
mkdir $BSRCL_DIR/round_$ROUND
cd $BSRCL_DIR/round_$ROUND

ls $BSRCL_DIR/round_2/*.vcf >interval_vcf.list

MERGE_JOB_NAME=$SAMPLE".merge_cohort_round_"$ROUND
THREADS=1

IN_LIST="$INTERVAL.intList"
OUT_VCF="round_$ROUND/$INTERVAL.vcf"
    
MERGE="gatk --java-options "'"-Xmx4G"'" MergeVcfs -I  $IN_LIST -O $OUT_VCF -R $REFERENCE"

MERGE_QSUB="$QSUB -pe mpi $THREADS -N $MERGE_JOB_NAME -o logs/$MERGE_JOB_NAME.log"
echo $MERGE >scripts/$MERGE_JOB_NAME.sh

cat scripts/$MERGE_JOB_NAME.sh | $MERGE_QSUB
 
#<<<<<<<<<<<<<<<<<<< wait till all finished >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# SORT COHORT VCF --------------------------------------------------------------
SORT_JOB_NAME=$SAMPLE".sort_indiv"
THREADS=12

IN_GVCF=$OUT_GVCF
OUT_GVCF="$SAMPLE.g.vcf"
    
SORT="$SINGULARITY gatk --java-options "'"-Xmx4G"'" SortVcf -I $IN_GVCF -O $OUT_GVCF"

SORT_QSUB="$QSUB -pe mpi $THREADS -N $SORT_JOB_NAME -o logs/$SORT_JOB_NAME.log -hold_jid $MERGE_JOB_NAME"
echo $SORT >scripts/$SORT_JOB_NAME.sh

cat scripts/$SORT_JOB_NAME.sh | $SORT_QSUB



#if i want to check for completion
#for LOG in *.o; do
#  if [[ $(grep "picard.vcf.MergeVcfs done" $LOG) ]]; then
#    PASSED=$((PASSED+1))
#  else
#    FAILED=$((FAILED+1))
#    echo $LOG
#    fi 
#    TOTAL=$((TOTAL+1))
#done

# CLEANUP ----------------------------------------------------------------------

# remove all unnecessary files
# put sorted cohort vcf somehwere safe
