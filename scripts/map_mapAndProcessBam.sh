#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# map_process_bam.sh - maps reads to genome and post-processes bam file for gatk

# TODO(nplatt): add code to clean up dir of tmp bam files
# TODO(nplatt): update comments

source /master/nplatt/sH_hybridization/scripts/set-env.sh

mkdir $MAP_DIR $MAP_DIR/logs $MAP_DIR/scripts

cd $MAP_DIR

#Map R1, R2, RX
for SAMPLE in $(cat $SAMPLE_LIST); do
    for READ in R1 R2 RX; do

        BWA_JOB_NAME=$SAMPLE"_"$READ".bwa_aln"

        IN_READS=$FILTER_READS_DIR/$SAMPLE"_filtered_"$READ".fastq.gz"
        OUT_SAI=$MAP_DIR/$SAMPLE"_"$READ".sai"
 
        THREADS=12

        BWA_QSUB="$QSUB -pe mpi $THREADS -N $BWA_JOB_NAME -o logs/$BWA_JOB_NAME.log -hold_jid bwa_index,$SAMPLE.filter_reads"
   
        BWA="$SINGULARITY bwa aln -t $THREADS -f $OUT_SAI $REFERENCE $IN_READS"
    
        echo $BWA >scripts/$BWA_JOB_NAME.sh
        cat scripts/$BWA_JOB_NAME.sh | $BWA_QSUB

    done
    
    # SAMPE --------------------------------------------------------------------
    SAMPE_JOB_NAME=$SAMPLE".sampe"
    THREADS=12

    IN_SAI_R1=$MAP_DIR/$SAMPLE"_R1.sai"
    IN_SAI_R2=$MAP_DIR/$SAMPLE"_R2.sai"
    IN_FQ_R1=$FILTER_READS_DIR/$SAMPLE"_filtered_R1.fastq.gz"
    IN_FQ_R2=$FILTER_READS_DIR/$SAMPLE"_filtered_R2.fastq.gz"
    
    OUT_BAM=$MAP_DIR/$SAMPLE"_samPE.bam"
    
    SAMPE="$SINGULARITY bwa sampe $REFERENCE $IN_SAI_R1 $IN_SAI_R2 $IN_FQ_R1 $IN_FQ_R2 | samtools view -Sb -F 4 - >$OUT_BAM"

    SAMPE_QSUB="$QSUB -pe mpi $THREADS -N $SAMPE_JOB_NAME -o logs/$SAMPE_JOB_NAME.log -hold_jid $SAMPLE"_R1.bwa_aln",$SAMPLE"_R2.bwa_aln
    echo $SAMPE >scripts/$SAMPE_JOB_NAME.sh



    # SORT_SAMPE ----------------------------------------------------------------
    SORTSAMPE_JOB_NAME=$SAMPLE".sort_sampe"
    THREADS=12

    IN_BAM=$MAP_DIR/$SAMPLE"_samPE.bam"
   
    OUT_PREFIX=$MAP_DIR/$SAMPLE"_samPE_sorted"
    
    SORTSAMPE="$SINGULARITY samtools sort $IN_BAM $OUT_PREFIX"

    SORTSAMPE_QSUB="$QSUB -pe mpi $THREADS -N $SORTSAMPE_JOB_NAME -o logs/$SORTSAMPE_JOB_NAME.log -hold_jid $SAMPE_JOB_NAME"
    echo $SORTSAMPE >scripts/$SORTSAMPE_JOB_NAME.sh

    # SAMSE --------------------------------------------------------------------
    SAMSE_JOB_NAME=$SAMPLE".samse"
    THREADS=12

    IN_SAI_RX=$MAP_DIR/$SAMPLE"_RX.sai"
    IN_FQ_RX=$FILTER_READS_DIR/$SAMPLE"_filtered_RX.fastq.gz"
    OUT_BAM=$MAP_DIR/$SAMPLE"_samSE.bam"

    SAMSE="$SINGULARITY bwa samse $REFERENCE $IN_SAI_RX $IN_FQ_RX | samtools view -Sb -F 4 - >$OUT_BAM"      

    SAMSE_QSUB="$QSUB -pe mpi $THREADS -N $SAMSE_JOB_NAME -o logs/$SAMSE_JOB_NAME.log -hold_jid $SAMPLE"_RX.bwa_aln
    echo $SAMSE >scripts/$SAMSE_JOB_NAME.sh


    # SORT_SAMSE ----------------------------------------------------------------
    SORTSAMSE_JOB_NAME=$SAMPLE".sort_samse"
    THREADS=12

    IN_BAM=$MAP_DIR/$SAMPLE"_samSE.bam"
   
    OUT_PREFIX=$MAP_DIR/$SAMPLE"_samSE_sorted"
    
    SORTSAMSE="$SINGULARITY samtools sort $IN_BAM $OUT_PREFIX"

    SORTSAMSE_QSUB="$QSUB -pe mpi $THREADS -N $SORTSAMSE_JOB_NAME -o logs/$SORTSAMSE_JOB_NAME.log -hold_jid $SAMSE_JOB_NAME"
    echo $SORTSAMSE >scripts/$SORTSAMSE_JOB_NAME.sh



    # MERGE --------------------------------------------------------------------
    MERGE_JOB_NAME="$SAMPLE.merge"
    THREADS=12

    IN_SAMSE=$MAP_DIR/$SAMPLE"_samSE_sorted.bam"
    IN_SAMPE=$MAP_DIR/$SAMPLE"_samPE_sorted.bam"

    OUT_MERGEDBAM=$MAP_DIR/$SAMPLE"_merged.bam"

    MERGE="$SINGULARITY samtools merge $OUT_MERGEDBAM $IN_SAMSE $IN_SAMPE"

    MERGE_QSUB="$QSUB -pe mpi $THREADS -N $MERGE_JOB_NAME -o logs/$MERGE_JOB_NAME.log -hold_jid $SORTSAMSE_JOB_NAME,$SORTSAMPE_JOB_NAME"
    echo $MERGE >scripts/$MERGE_JOB_NAME.sh



    # RGS ----------------------------------------------------------------------
    RGS_JOB_NAME="$SAMPLE.add_readgroups"
    THREADS=1

    IN_FASTQ=$RAW_READS_DIR/$SAMPLE/$SAMPLE"_R1.fastq.gz"
    IN_BAM=$MAP_DIR/$SAMPLE"_merged.bam"
    OUT_BAM=$MAP_DIR/$SAMPLE"_merged_RGs.bam"

    HEADER_LINE=$(zcat $IN_FASTQ | head -n 1)

    CELL=$(echo $HEADER_LINE | cut -f3 -d":")
    LANE=$(echo $HEADER_LINE | cut -f4 -d":")
    INDEX=$(echo $HEADER_LINE | cut -f10 -d":")

    RGS="$SINGULARITY gatk AddOrReplaceReadGroups \
        --INPUT=$IN_BAM \
        --OUTPUT=$OUT_BAM \
        --RGID=$CELL.$LANE \
        --RGLB=library1 \
        --RGPL=illumina \
        --RGPU=$CELL.$INDEX.$LANE \
        --RGSM=$SAMPLE"

    RGS_QSUB="$QSUB -pe mpi $THREADS -N $RGS_JOB_NAME -o logs/$RGS_JOB_NAME.log -hold_jid $MERGE_JOB_NAME"
    echo $RGS >scripts/$RGS_JOB_NAME.sh



    # SORT- --------------------------------------------------------------------
    SORT_JOB_NAME="$SAMPLE.sort"
    THREADS=12

    IN_BAM=$MAP_DIR/$SAMPLE"_merged_RGs.bam"

    OUT_BAM=$MAP_DIR/$SAMPLE"_merged_RGs_sorted.bam"

    SORT="$SINGULARITY gatk SortSam --INPUT $IN_BAM --OUTPUT $OUT_BAM --SORT_ORDER=coordinate"
    
    SORT_QSUB="$QSUB -pe mpi $THREADS -N $SORT_JOB_NAME -o logs/$SORT_JOB_NAME.log -hold_jid $RGS_JOB_NAME"
    echo $SORT >scripts/$SORT_JOB_NAME.sh

 

    # DUPES --------------------------------------------------------------------
    DUPES_JOB_NAME="$SAMPLE.mark_duplicates"
    THREADS=12

    IN_BAM=$MAP_DIR/$SAMPLE"_merged_RGs_sorted.bam"

    OUT_BAM=$MAP_DIR/$SAMPLE"_processed.bam"
    OUT_METRICS=$MAP_DIR/$SAMPLE"_nodupes.metrics"

    DUPES="$SINGULARITY gatk MarkDuplicates \
            --INPUT $IN_BAM \
            --OUTPUT $OUT_BAM \
            --METRICS_FILE $OUT_METRICS \
            --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 900"
    
    DUPES_QSUB="$QSUB -pe mpi $THREADS -N $DUPES_JOB_NAME -o logs/$DUPES_JOB_NAME.log -hold_jid $SORT_JOB_NAME"
    echo $DUPES >scripts/$DUPES_JOB_NAME.sh
    

    # INDEX --------------------------------------------------------------------
    INDEX_JOB_NAME="$SAMPLE.index"
    THREADS=12

    IN_BAM=$MAP_DIR/$SAMPLE"_processed.bam"

    INDEX="$SINGULARITY samtools index $IN_BAM"
    
    INDEX_QSUB="$QSUB -pe mpi $THREADS -N $INDEX_JOB_NAME -o logs/$INDEX_JOB_NAME.log -hold_jid $DUPES_JOB_NAME"
    echo $INDEX >scripts/$INDEX_JOB_NAME.sh


    # CLEANUP-------------------------------------------------------------------
    #CLEANUP_JOB_NAME="$SAMPLE.clean_up"
    #THREADS=1

    #CLEANUP="rm $MAP_DIR/$SAMPLE"_R1.sai "$MAP_DIR/$SAMPLE"_R2.sai "$MAP_DIR/$SAMPLE"_RX.sai "$MAP_DIR/$SAMPLE"_samSE.bam "$MAP_DIR/$SAMPLE"_samPE.bam "$MAP_DIR/$SAMPLE"_samSE_sorted.bam "$MAP_DIR/$SAMPLE"_samPE_sorted.bam "$MAP_DIR/$SAMPLE"_merged.bam "$MAP_DIR/$SAMPLE"_merged_RGs.bam "$MAP_DIR/$SAMPLE"_merged_RGs_sorted.bam; chmod a-w "$MAP_DIR/$SAMPLE"_processed.bam"
    #    "

    #CLEANUP_QSUB="$QSUB -pe mpi $THREADS -N $CLEANUP_JOB_NAME -o logs/$CLEANUP_JOB_NAME.log -hold_jid $DUPES_JOB_NAME"
    #echo $CLEANUP >scripts/$CLEANUP_JOB_NAME.sh



    # SUBMIT JOBS --------------------------------------------------------------
    cat scripts/$SAMPE_JOB_NAME.sh     | $SAMPE_QSUB
    cat scripts/$SORTSAMPE_JOB_NAME.sh | $SORTSAMPE_QSUB
    cat scripts/$SAMSE_JOB_NAME.sh     | $SAMSE_QSUB
    cat scripts/$SORTSAMSE_JOB_NAME.sh | $SORTSAMSE_QSUB
    cat scripts/$MERGE_JOB_NAME.sh     | $MERGE_QSUB
    cat scripts/$RGS_JOB_NAME.sh       | $RGS_QSUB
    cat scripts/$SORT_JOB_NAME.sh      | $SORT_QSUB
    cat scripts/$DUPES_JOB_NAME.sh     | $DUPES_QSUB
    cat scripts/$INDEX_JOB_NAME.sh | $INDEX_QSUB
    #cat scripts/$CLEANUP_JOB_NAME.sh | $CLEANUP_QSUB

done



