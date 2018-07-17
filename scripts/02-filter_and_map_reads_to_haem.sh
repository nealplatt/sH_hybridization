#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh

##DOWNLOAD SRA DATA
#cd $SEQ_DIR
#for SRA_ACCESSION in ${SRA_SAMPLES[@]}; do
#    ${ENVIRONMENTS["SINGULARITY"]} fastq-dump --split-files --gzip $SRA_ACCESSION &
#done

#wait

#rename _1.fastq.gz _R1.fastq.gz *_1.fastq.gz
#rename _2.fastq.gz _R2.fastq.gz *_2.fastq.gz

#START FILTERING/MAPPING PROCESS
cd $RESULTS_DIR
for SAMPLE in "${SAMPLES[@]}"; do

    echo"" >$SAMPLE"_process.sh"

    ############ FILTER_READS
    JID=$SAMPLE"_filter"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"     
    HOLD="-hold_jid "$SAMPLE"_R1_map,"$SAMPLE"_R2_map"

    R1=$SEQ_DIR/$SAMPLE"_R1.fastq.gz"
    R2=$SEQ_DIR/$SAMPLE"_R1.fastq.gz"
    PE_R1=$FILTER_DIR/$SAMPLE"_filtered_paired_R1.fastq.gz"
    PE_R2=$FILTER_DIR/$SAMPLE"_filtered_paired_R2.fastq.gz"
    SE_R1=$FILTER_DIR/$SAMPLE"_filtered_unpaired_R1.fastq.gz"
    SE_R2=$FILTER_DIR/$SAMPLE"_filtered_unpaired_R2.fastq.gz"
    SE=$FILTER_DIR/$SAMPLE"_filtered_unpaired_RX.fastq.gz"
    
    CMD="${ENVIRONMENTS[$ENV]} trimmomatic \
            PE \
            -threads $THREADS \
            -phred33 \
            $R1 $R2 \
            $PE_R1 $SE_R1 $PE_R2 $SE_R2 \
            LEADING:10 \
            TRAILING:10 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36;
          
        zcat $SE_R1 $SE_R2 | gzip >$SE;
        rm $SE_R1 $SE_R2"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"
      
    ############ MAP R1 and R2 READS
    for READ in R1 R2; do
        JID=$SAMPLE"_"$READ"_map"
        LOG=$LOGS_DIR/$JID".log"
        THREADS=12 
        ENV="CONDA_SNPS"
        HOLD="-hold_jid "$SAMPLE"_filter"
    
        SAI=$MAP_DIR/$SAMPLE"_"$READ".sai"
        FQ=$FILTER_DIR/$SAMPLE"_filtered_paired_"$READ".fastq.gz"
        GENOME=$HAE_GENOME

        JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

        CMD="${ENVIRONMENTS[$ENV]} bwa aln -t 12 -f $SAI $GENOME $FQ"
            
        echo $CMD | $JOB_QSUB
        echo $CMD >>$SAMPLE"_process.sh"
    done

    ############ MAP RX READS
    JID=$SAMPLE"_RX_map"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12 
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_filter"
    
    SAI=$MAP_DIR/$SAMPLE"_RX.sai"
    FQ=$FILTER_DIR/$SAMPLE"_filtered_unpaired_RX.fastq.gz"
    GENOME=$HAE_GENOME

    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    CMD="${ENVIRONMENTS[$ENV]} bwa aln -t 12 -f $SAI $GENOME $FQ"
            
    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############ SAMPE
    JID=$SAMPLE"_sampe"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"     
    HOLD="-hold_jid "$SAMPLE"_R1_map,"$SAMPLE"_R2_map"

    GENOME=$HAE_GENOME
    SAI_1=$MAP_DIR/$SAMPLE"_R1.sai"
    SAI_2=$MAP_DIR/$SAMPLE"_R2.sai"
    PE_1=$FILTER_DIR/$SAMPLE"_filtered_paired_R1.fastq.gz"
    PE_2=$FILTER_DIR/$SAMPLE"_filtered_paired_R2.fastq.gz"
    SAMPE_BAM=$MAP_DIR/$SAMPLE"_sampe.bam"
    
    CMD="${ENVIRONMENTS[$ENV]} bwa sampe $GENOME $SAI_1 $SAI_2 $PE_1 $PE_2 | samtools view -Sb -F 4 - >$SAMPE_BAM"      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############ SAMSE
    JID=$SAMPLE"_samse"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_RX_map"

    GENOME=$HAE_GENOME
    SAI_X=$MAP_DIR/$SAMPLE"_RX.sai"
    PE_X=$FILTER_DIR/$SAMPLE"_filtered_unpaired_RX.fastq.gz"
    SAMSE_BAM=$MAP_DIR/$SAMPLE"_samse.bam"
    
    CMD="${ENVIRONMENTS[$ENV]} bwa samse $GENOME $SAI_X $PE_X | samtools view -Sb -F 4 - >$SAMSE_BAM"      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############ SORT SAMPE
    JID=$SAMPLE"_sampe_sort"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_sampe"

    SAMPE_BAM=$MAP_DIR/$SAMPLE"_sampe.bam"
    SORTED_SAMPE=$MAP_DIR/$SAMPLE"_sorted_sampe.bam"
    
    CMD="${ENVIRONMENTS[$ENV]}  samtools sort --threads $THREADS -o $SORTED_SAMPE $SAMPE_BAM"      
    JOB_QSUB="$QSUB -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############ SORT SAMSE
    JID=$SAMPLE"_samse_sort"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_samse"

    SAMSE_BAM=$MAP_DIR/$SAMPLE"_samse.bam"
    SORTED_SAMSE=$MAP_DIR/$SAMPLE"_sorted_samse.bam"
    
    CMD="${ENVIRONMENTS[$ENV]} samtools sort --threads $THREADS -o $SORTED_SAMSE $SAMSE_BAM "      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############ MERGE SAMSE AND SAMSE
    JID=$SAMPLE"_merge"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_samse_sort,"$SAMPLE"_sampe_sort"

    SORTED_SAMSE=$MAP_DIR/$SAMPLE"_sorted_samse.bam"
    SORTED_SAMPE=$MAP_DIR/$SAMPLE"_sorted_sampe.bam"
    MERGED_BAM=$MAP_DIR/$SAMPLE".bam"    

    CMD="${ENVIRONMENTS[$ENV]} samtools merge $MERGED_BAM $SORTED_SAMSE $SORTED_SAMPE"      
    
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############ ADD READ GROUP INFO
    JID=$SAMPLE"_add_rgs"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="SINGULARITY"
    HOLD="-hold_jid "$SAMPLE"_merge"
    
    MERGED_BAM=$MAP_DIR/$SAMPLE".bam"
    RG_BAM=$MAP_DIR/$SAMPLE"_rg.bam"        

    LANE=${RG_LANE[$SAMPLE]}
    CELL=${RG_CELL[$SAMPLE]}
    INDEX=${RG_INDEX[$SAMPLE]}

    CMD="${ENVIRONMENTS[$ENV]} \
            gatk AddOrReplaceReadGroups \
                --INPUT=$MERGED_BAM \
                --OUTPUT=$RG_BAM \
                --RGID=$CELL.$LANE \
                --RGLB=library1 \
                --RGPL=illumina \
                --RGPU=$CELL.$INDEX.$LANE  \
                --RGSM=$SAMPLE"
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############ MARK DUPLICATE READS
    JID=$SAMPLE"_mark_dups"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="SINGULARITY"
    HOLD="-hold_jid "$SAMPLE"_add_rgs"
    
    RG_BAM=$MAP_DIR/$SAMPLE"_rg.bam"        
    FINAL_BAM=$MAP_DIR/$SAMPLE"_processed.bam"        
    METRICS=$MAP_DIR/$SAMPLE"_processed.metrics"        

    CMD="${ENVIRONMENTS[$ENV]} \
            gatk MarkDuplicates \
                --INPUT $RG_BAM \
                --OUTPUT $FINAL_BAM \
                --METRICS_FILE $METRICS \
                --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 900"
 
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############ INDEX BAM FILE
    JID=$SAMPLE"_bam_index"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=1
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_mark_dups"
    
    FINAL_BAM=$MAP_DIR/$SAMPLE"_processed.bam"        

    CMD="${ENVIRONMENTS[$ENV]} samtools index $FINAL_BAM"
 
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"
   
done

#once all jobs are done, clean up dir
for SAMPLE in "${SAMPLES[@]}"; do
    rm $SAMPLE".bam"
    rm $SAMPLE"_R1.sai"
    rm $SAMPLE"_R2.sai"
    rm $SAMPLE"_rg.bam"
    rm $SAMPLE"_RX.sai"
    rm $SAMPLE"_sampe.bam"
    rm $SAMPLE"_samse.bam"
    rm $SAMPLE"_sorted_sampe.bam"
    rm $SAMPLE"_sorted_samse.bam"
    rm $SAMPLE"__noDupes.metrics"
done

mkdir metrics

