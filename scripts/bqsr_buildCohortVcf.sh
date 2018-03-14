#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr_buildCohortVcf.sh - after genotyping, build and filter cohort vcf for 
#       recalibration

# TODO(nplatt): update comments

WAIT_FOR_CLEAR_QUEUE

cd $ROUND"_cohort_vcf"

# MERGE ROUND 1 ----------------------------------------------------------------
ls $BQSR_DIR/$ROUND"_genotype"/*.vcf >interval_vcf.list

$SINGULARITY split -d -n l/1000 --additional-suffix .list interval_vcf.list m1.

for INTERVAL in $(seq -w 0 999); do

    JOB_NAME="snp.merge_cohort_"$INTERVAL"_"$ROUND
    THREADS=1
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND=""
    SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

    IN_LIST="m1."$INTERVAL".list"
    OUT_VCF="m1.$INTERVAL.vcf"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    CMD="$SINGULARITY gatk MergeVcfs -I $IN_LIST -O $OUT_VCF -R $REFERENCE"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"

done


#
#                               <...wait...>
#
#sleep while all jobs are running
echo "waiting for queue to clear"
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
        JOB_NAME="snp.merge_cohort_"$INTERVAL"_"$ROUND
        THREADS=1
        LOG="$LOGS_DIR/$JOB_NAME.log" 
        DEPEND=""
        SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"
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
echo "waiting for queue to clear"
WAIT_FOR_CLEAR_QUEUE
#
#                               <...wait...>
#

# MERGE ROUND 2 ----------------------------------------------------------------
rm *.list
ls m1.*.vcf >interval_vcf.list

$SINGULARITY split -d -n l/100 --additional-suffix .list interval_vcf.list m2.

for INTERVAL in $(seq -w 0 99); do

    JOB_NAME="snp.merge2_cohort_"$INTERVAL"_"$ROUND
    THREADS=1
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND=""
    SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

    IN_LIST="$WORK_DIR/m2."$INTERVAL".list"
    OUT_VCF="$WORK_DIR/m2.$INTERVAL.vcf"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    CMD="$SINGULARITY gatk MergeVcfs -I $IN_LIST -O $OUT_VCF -R $REFERENCE"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"

done


#
#                               <...wait...>
#
#sleep while all jobs are running
echo "waiting for queue to clear"
WAIT_FOR_CLEAR_QUEUE
#
#                               <...wait...>
#


# MERGE ROUND 3 ----------------------------------------------------------------
rm *.list
ls m2.*.list >interval_vcf.list

JOB_NAME="snp.merge3_cohort_"$ROUND
THREADS=1
LOG="$LOGS_DIR/$JOB_NAME.log" 
DEPEND=""
SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

IN_LIST="interval_vcf.list"
OUT_VCF="tmp_cohort.vcf"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

CMD="$SINGULARITY gatk MergeVcfs -I $IN_LIST -O $OUT_VCF -R $REFERENCE"

DELETE $LOG $SCRIPT
SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"

# SORT COHORT VCF --------------------------------------------------------------
SORT_JOB_NAME="snp.sort_cohort"
THREADS=12
LOG="$LOGS_DIR/$JOB_NAME.log"
DEPEND="-hold_jid snp.merge3_cohort_"$ROUND
SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

IN_GVCF="tmp_cohort.vcf"
OUT_GVCF=$ROUND"_vcfs/cohort_raw_"$ROUND.vcf

CMD="$SINGULARITY gatk --java-options "'"-Xmx8G"'" SortVcf -I $IN_GVCF -O $OUT_GVCF"

SORT_QSUB="$QSUB -pe mpi $THREADS -N $SORT_JOB_NAME -o logs/$SORT_JOB_NAME.log"

DELETE $LOG $SCRIPT
SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"

# SELECT SNPS ------------------------------------------------------------------
JOB_NAME="snp.cohort_select_snps"
THREADS=1
LOG="$LOGS_DIR/$JOB_NAME.log" 
DEPEND="snp.sort_cohort"
SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

IN_VCF=$ROUND"_vcfs/cohort_raw_"$ROUND.vcf
OUT_VCF=$ROUND"_vcfs/cohort_"$ROUND"_rawSNPS.g.vcf"
    
CMD="$SINGULARITY gatk SelectVariants \
    -V $IN_VCF \
    -select-type SNP \
    -O $OUT_VCF \
    -R $REFERENCE"

DELETE $LOG $SCRIPT
SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"

# FILTER SNPS ------------------------------------------------------------------
JOB_NAME="snp.cohort_filter_snps"
THREADS=1
LOG="$LOGS_DIR/$JOB_NAME.log" 
DEPEND="-hold_jid snp.cohort_select_snps"
SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

IN_VCF=$ROUND"_vcfs/cohort_"$ROUND"_rawSNPS.g.vcf"
OUT_VCF=$ROUND"_vcfs/cohort_"$ROUND"_filteredSNPS.g.vcf"

CMD="$SINGULARITY gatk VariantFiltration \
    -R $REFERENCE \
    -V $IN_VCF \
    --filter-expression "'"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"'" \
    --filter-name "'"recommended_snp_filter"'" \
    -O $OUT_VCF"

DELETE $LOG $SCRIPT
SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"

# SELECT INDELS ----------------------------------------------------------------
JOB_NAME="snp.cohort_select_indels"
THREADS=1
LOG="$LOGS_DIR/$JOB_NAME.log" 
DEPEND="snp.sort_cohort"
SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

IN_VCF="$WORK_DIR/cohort_$ROUND.g.vcf"
OUT_VCF=$ROUND"_vcfs/cohort_"$ROUND"_rawINDELS.g.vcf"
    
CMD="$SINGULARITY gatk SelectVariants \
    -V $IN_VCF \
    -select-type INDEL \
    -O $OUT_VCF \
    -R $REFERENCE"

DELETE $LOG $SCRIPT
SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"


# FILTER INDELS ------------------------------------------------------------------
JOB_NAME="snp.cohort_filter_indels"
THREADS=1
LOG="$LOGS_DIR/$JOB_NAME.log" 
DEPEND="-hold_jid snp.cohort_select_indels"
SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

IN_VCF=$ROUND"_vcfs/cohort_"$ROUND"_rawINDELS.g.vcf"
OUT_VCF=$ROUND"_vcfs/cohort_"$ROUND"_filteredINDELS.g.vcf"

CMD="$SINGULARITY gatk VariantFiltration \
    -R $REFERENCE \
    -V $IN_VCF \
    --filter-expression "'"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"'" \
    --filter-name "'"recommended_indel_filter"'" \
    -O $OUT_VCF"

DELETE $LOG $SCRIPT
SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"


# MERGE VARIANTS----------------------------------------------------------------
JOB_NAME="snp.merge_filtered_variants"
THREADS=1
LOG="$LOGS_DIR/$JOB_NAME.log" 
DEPEND="-hold_jid snp.cohort_filter_indels,snp.cohort_filter_snps"
SCRIPT="$SCRIPTS_DIR/$JOB_NAME.sh"

JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

IN_SNP_VCF=$ROUND"_vcfs/cohort_"$ROUND"_filteredSNPS.g.vcf"
IN_INDEL_VCF=$ROUND"_vcfs/cohort_"$ROUND"_filteredINDELS.g.vcf"
OUT_VCF=$ROUND"_vcfs/cohort_"$ROUND"_filteredVariants.g.vcf"

CMD="$SINGULARITY gatk CombineVariants \
   -R $REFERENCE \
   --variant:snps $IN_SNP_VCF \
   --variant:indels $IN_INDEL_VCF \
   -O $OUT_VCF \
   -genotypeMergeOptions PRIORITIZE \
   -priority snps,indels"

DELETE $LOG $SCRIPT
SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"


################################################################################
#                       R E C A L I B R A T I O N 
################################################################################

