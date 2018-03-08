#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# build_gatk_intervals.sh - identifies (indirectly) regions in the Sh. genome
#  with no mapped reads, to split genome for downstream analyses

# TODO(nplatt): Update comments

#load variables
source master/nplatt/sH_hybridizationscripts/set_env.sh

mkdir $INTERVALS_DIR $INTERVALS_DIR/scripts $INTERVALS_DIR/logs

cd $INTERVALS_DIR

#-------------------------------------------------------------------------------
# Create a list of mapping intervals rather than whole contigs

for SAMPLE in $(cat $SAMPLE_LIST); do
    # BED --------------------------------------------------------------------
    BED_JOB_NAME=$SAMPLE".bed"
    THREADS=4

    IN_BAM=$MAP_DIR/$SAMPLE"_processed.bam"
    
    OUT_BED="tmp.$SAMPLE.bed"
    
    BED="$SINGULARITY bedtools bamtobed -i $IN_BAM >$OUT_BED"

    BED_QSUB="$QSUB -pe mpi $THREADS -N $BED_JOB_NAME -o logs/$BED_JOB_NAME.log"
    echo $BED >scripts/$BED_JOB_NAME.sh


    # SORT ---------------------------------------------------------------------
    SORT_JOB_NAME=$SAMPLE".sort"
    THREADS=12

    IN_BED=$OUT_BED
    
    OUT_BED="tmp.$SAMPLE.sorted.bed"
    
    SORT="$SINGULARITY bedtools sort -i $IN_BED >$OUT_BED" 

    SORT_QSUB="$QSUB -pe mpi $THREADS -N $SORT_JOB_NAME -o logs/$sort_JOB_NAME.log -hold_jid $BED_JOB_NAME"
    echo $SORT >scripts/$SORT_JOB_NAME.sh


    # MERGE --------------------------------------------------------------------
    MERGE_JOB_NAME=$SAMPLE".merge"
    THREADS=4

    IN_BED=$OUT_BED
    
    OUT_BED="tmp.$SAMPLE.merged.bed"
    
    MERGE="$SINGULARITY bedtools merge -d 0 -i $IN_BED >$OUT_BED" 

    MERGE_QSUB="$QSUB -pe mpi $THREADS -N $MERGE_JOB_NAME -o logs/$MERGE_JOB_NAME.log -hold_jid $SORT_JOB_NAME"
    echo $MERGE >scripts/$MERGE_JOB_NAME.sh
    

    # submit--------------------------------------------------------------------
    cat scripts/$BED_JOB_NAME.sh | $BED_QSUB
    cat scripts/$SORT_JOB_NAME.sh | $SORT_QSUB
    cat scripts/$MERGE_JOB_NAME.sh | $MERGE_QSUB
  

done

#CHECK FOR COMPLETION
for SAMPLE in $(cat $SAMPLE_LIST); do
    #Check to see if these files exist
    for FILE in "$SAMPLE.bed" "tmp.$SAMPLE.sorted.bed" "tmp.$SAMPLE.merged.bed"; do
        if [ ! -f $FILE ]; then
            echo "$SAMPLE"
        fi
    done
done

#any samples with missing files need to be resubmitted





#combine all bed files when done
# CATALL------------------------------------------------------------------------
CATALL_JOB_NAME="catall"
THREADS=1

OUT_BED="tmp.cat.bed"

CATALL="cat $INTERVALS_DIR/tmp.*.merged.bed >$OUT_BED"

CATALL_QSUB="$QSUB -pe mpi $THREADS -N $CATALL_JOB_NAME -o logs/$CATALL_JOB_NAME.log"
echo $CATALL >scripts/$CATALL_JOB_NAME.sh


# SORTALL-----------------------------------------------------------------------
SORTALL_JOB_NAME="sortall"
THREADS=12

IN_BED=$OUT_BED
OUT_BED="tmp.sort.bed"

SORTALL="$SINGULARITY bedtools sort -i $IN_BED >$OUT_BED"

SORTALL_QSUB="$QSUB -pe mpi $THREADS -N $SORTALL_JOB_NAME -o logs/$SORTALL_JOB_NAME.log -hold_jid $CATALL_JOB_NAME"
echo $SORTALL >scripts/$SORTALL_JOB_NAME.sh


# SLOPALL-----------------------------------------------------------------------
SLOPALL_JOB_NAME="slopall"
THREADS=12

IN_BED=$OUT_BED
OUT_BED="tmp.slop.bed"

SLOPALL="$SINGULARITY bedtools slop -b 150 -g $REFERENCE.fai -i $IN_BED >$OUT_BED"

SLOPALL_QSUB="$QSUB -pe mpi $THREADS -N $SLOPALL_JOB_NAME -o logs/$SLOPALL_JOB_NAME.log -hold_jid $SORTALL_JOB_NAME"
echo $SLOPALL >scripts/$SLOPALL_JOB_NAME.sh


# MERGEALL ---------------------------------------------------------------------
MERGEALL_JOB_NAME="mergeall"
THREADS=1

IN_BED=$OUT_BED
OUT_BED="all_intervals.bed"

MERGEALL="$SINGULARITY bedtools merge -d 200 -c 1 -o count -i $IN_BED >$OUT_BED"

MERGEALL_QSUB="$QSUB -pe mpi $THREADS -N $MERGEALL_JOB_NAME -o logs/$MERGEALL_JOB_NAME.log -hold_jid $SLOPALL_JOB_NAME"
echo $MERGEALL >scripts/$MERGEALL_JOB_NAME.sh


# FILTALL ----------------------------------------------------------------------
FILTALL_JOB_NAME="filtall"
THREADS=1

IN_BED=$OUT_BED
OUT_LIST="$INTERVALS_DIR/all_filtered_intervals.bed"

FILTALL="$SINGULARITY awk '"'$4'" >= 24 {print "'$1+1":"$2"-"$3'"}' $IN_BED >$OUT_LIST;"

FILTALL_QSUB="$QSUB -pe mpi $THREADS -N $FILTALL_JOB_NAME -o logs/$FILTALL_JOB_NAME.log -hold_jid $MERGEALL_JOB_NAME"
echo $FILTALL >scripts/$FILTALL_JOB_NAME.sh


# SPLITALL ---------------------------------------------------------------------
SPLITALL_JOB_NAME="splitall"
THREADS=1

IN_LIST=$OUT_LIST

SPLITALL="$SINGULARITY split -n l/50 --numeric-suffixes --additional-suffix .list $IN_LIST filtered_interval.part"

SPLITALL_QSUB="$QSUB -pe mpi $THREADS -N $SPLITALL_JOB_NAME -o logs/$SPLITALL_JOB_NAME.log -hold_jid $FILTALL_JOB_NAME"
echo $SPLITALL >scripts/$SPLITALL_JOB_NAME.sh


# submit--------------------------------------------------------------------
cat scripts/$CATALL_JOB_NAME.sh   | $CATALL_QSUB
cat scripts/$SORTALL_JOB_NAME.sh  | $SORTALL_QSUB
cat scripts/$SLOPALL_JOB_NAME.sh  | $SLOPALL_QSUB
cat scripts/$MERGEALL_JOB_NAME.sh | $MERGEALL_QSUB
cat scripts/$FILTALL_JOB_NAME.sh  | $FILTALL_QSUB
cat scripts/$SPLITALL_JOB_NAME.sh | $SPLITALL_QSUB  



# CREATE CONTIG LIST JIC -------------------------------------------------------
grep ">" $REFERENCE | cut -f1 -d" " | sed 's/>//' >contigs.list


# CLEAN ------------------------------------------------------------------------
rm tmp.*




