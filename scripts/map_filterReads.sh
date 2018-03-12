#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# filter_reads.sh - filter all reads with trimmomatic

# TODO(nplatt): update comments

#load variables
source /master/nplatt/sH_hybridization/scripts/set-env.sh

#Filter reads wit trimmomatic
mkdir $FILTER_READS_DIR $FILTER_READS_DIR/scripts $FILTER_READS_DIR/logs

cd $FILTER_READS_DIR

for SAMPLE in $(cat $SAMPLE_LIST); do

    JOB_NAME=$SAMPLE.filter_reads
   
    IN_R1=$RAW_READS_DIR/$SAMPLE/$SAMPLE"_R1.fastq.gz"
    IN_R2=$RAW_READS_DIR/$SAMPLE/$SAMPLE"_R2.fastq.gz"

    OUT_R1_PE=$FILTER_READS_DIR/$SAMPLE"_filtered_R1.fastq.gz"
    OUT_R2_PE=$FILTER_READS_DIR/$SAMPLE"_filtered_R2.fastq.gz"
    OUT_R1_SE=$FILTER_READS_DIR/$SAMPLE"_filtered_SE_R1.fastq.gz"
    OUT_R2_SE=$FILTER_READS_DIR/$SAMPLE"_filtered_SE_R2.fastq.gz"
    OUT_SE=$FILTER_READS_DIR/$SAMPLE"_filtered_RX.fastq.gz"

    THREADS=12

    FILTER_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o logs/$JOB_NAME.log"

    CMD="$SINGULARITY java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar \
        PE \
        -threads $THREADS \
        -phred33 \
        $IN_R1 \
        $IN_R2 \
        $OUT_R1_PE \
        $OUT_R1_SE \
        $OUT_R2_PE \
        $OUT_R2_SE \
        LEADING:10 \
        TRAILING:10 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36

        zcat $OUT_R1_SE $OUT_R2_SE | gzip >$OUT_SE

        rm $OUT_R1_SE $OUT_R2_SE
        "   

    echo "$CMD" >scripts/$JOB_NAME.sh

    cat scripts/$JOB_NAME.sh | $FILTER_QSUB

done




