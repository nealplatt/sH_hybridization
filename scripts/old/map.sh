#BOV  - ERR119622 ERR103048
#CUR  - ERR119623

QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y"    
WORK_DIR="/master/nplatt/schisto_hybridization"
LOGS_DIR=$WORK_DIR/results/logs


declare -A ENVIRONMENT
ENVIRONMENT=( ["SINGULARITY"]="singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img"
              ["CONDA"]="source activate snp_calling;")

declare -A RG_CELL
declare -A RG_LANE
declare -A RG_INDEX

RG_CELL=( ["ERR119622"]="4" ["ERR103048"]="3" ["ERR119623"]="5" ["Sm.BR_1278.1"]="C3E5HACXX" ["Sm.BR_0447.1"]="C3E5HACXX" ["Sm.BR_2039.1"]="C3E5HACXX")
RG_LANE=( ["ERR119622"]="3" ["ERR103048"]="4" ["ERR119623"]="2" ["Sm.BR_1278.1"]="2" ["Sm.BR_0447.1"]="1" ["Sm.BR_2039.1"]="1")
RG_INDEX=( ["ERR119622"]="ERR119622" ["ERR103048"]="ERR103048" ["ERR119623"]="ERR119623cow" ["Sm.BR_1278.1"]="CAATGGAA" ["Sm.BR_0447.1"]="CTGTAGCC" ["Sm.BR_2039.1"]="AACGCTTA")


declare -A GENOMES
BOVIS_GENOME=$WORK_DIR"/data/genome/schBov_v1.fa"
CURRASONI_GENOME=$WORK_DIR"/data/genome/schCur_v1.fa"
MANSONI_GENOME=$WORK_DIR"/data/genome/schMan_v7.fa"

GENOMES=( ["ERR119622"]=$BOVIS_GENOME
          ["ERR103048"]=$BOVIS_GENOME 
          ["ERR119623"]=$CURRASONI_GENOME
          ["Sm.BR_1278.1"]=$MANSONI_GENOME 
          ["Sm.BR_0447.1"]=$MANSONI_GENOME 
          ["Sm.BR_2039.1"]=$MANSONI_GENOME )

#mapping reads to bovis genome
mkdir $WORK_DIR/results/01-processed_reads/outgroup_mapped_reads/
for SAMPLE in ERR119622 ERR103048; do

    touch $SAMPLE"_process.sh"

    for READ in R1 R2 RX; do
        JID=$SAMPLE"_"$READ
        LOG=$LOGS_DIR/$JID".log"
        THREADS=12 
        ENV="CONDA"
    
        SAI=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_"$READ".sai"
        FQ=$WORK_DIR"/results/01-processed_reads/filtered_reads/"$SAMPLE"_filtered_"*"paired_"$READ".fastq.gz"
        GENOME=${GENOMES[$SAMPLE]}

        JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"

        CMD="${ENVIRONMENT[$ENV]} bwa aln -t 12 -f $SAI $GENOME $FQ"
            
        #echo $CMD | $JOB_QSUB
        echo $CMD >>$SAMPLE"_process.sh"
    done
done

#mapping reads to currasoni genome
for SAMPLE in ERR119623; do
    for READ in R1 R2 RX; do
        JID=$SAMPLE"_"$READ
        LOG=$LOGS_DIR/$JID".log"
        THREADS=12 
        ENV="CONDA"     

        SAI=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_"$READ".sai"
        FQ=$WORK_DIR"/results/01-processed_reads/filtered_reads/"$SAMPLE"_filtered_"*"paired_"$READ".fastq.gz"
        GENOME=${GENOMES[$SAMPLE]}


        JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"
        CMD="${ENVIRONMENT[$ENV]} bwa aln -t 12 -f $SAI $GENOME $FQ"
        
        #echo $CMD | $JOB_QSUB
        echo $CMD >>$SAMPLE"_process.sh"
    done
done

#mapping reads to mansoni genome
for SAMPLE in Sm.BR_1278.1 Sm.BR_0447.1 Sm.BR_2039.1; do
    for READ in R1 R2 RX; do
        JID=$SAMPLE"_"$READ
        LOG=$LOGS_DIR/$JID".log"
        THREADS=12 
        ENV="CONDA"     

        SAI=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_"$READ".sai"
        FQ=$WORK_DIR"/results/01-processed_reads/filtered_reads/"$SAMPLE"_filtered_"*"paired_"$READ".fastq.gz"
        GENOME=${GENOMES[$SAMPLE]}

        JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"
        CMD="${ENVIRONMENT[$ENV]} bwa aln -t 12 -f $SAI $GENOME $FQ" 

        #echo $CMD | $JOB_QSUB
        echo $CMD >>$SAMPLE"_process.sh"
    done
done

#processing mapped reads
for SAMPLE in ERR119622 ERR103048 ERR119623 Sm.BR_1278.1 Sm.BR_0447.1 Sm.BR_2039.1; do

    ############
    JID=$SAMPLE"_sampe"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA"     
    HOLD="-hold_jid "$SAMPLE"_R1,"$SAMPLE"_R2"

    GENOME=${GENOMES[$SAMPLE]}
    SAI_1=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_R1.sai"
    SAI_2=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_R2.sai"
    PE_1=$WORK_DIR"/results/01-processed_reads/filtered_reads/"$SAMPLE"_filtered_paired_R1.fastq.gz"
    PE_2=$WORK_DIR"/results/01-processed_reads/filtered_reads/"$SAMPLE"_filtered_paired_R2.fastq.gz"
    SAMPE_BAM=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_sampe.bam"
    
    CMD="${ENVIRONMENT[$ENV]} bwa sampe $GENOME $SAI_1 $SAI_2 $PE_1 $PE_2 | samtools view -Sb -F 4 - >$SAMPE_BAM"      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_samse"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA"
    HOLD="-hold_jid "$SAMPLE"_RX"

    GENOME=${GENOMES[$SAMPLE]}
    SAI_X=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_RX.sai"
    PE_X=$WORK_DIR"/results/01-processed_reads/filtered_reads/"$SAMPLE"_filtered_unpaired_RX.fastq.gz"
    SAMSE_BAM=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_samse.bam"
    
    CMD="${ENVIRONMENT[$ENV]} bwa samse $GENOME $SAI_X $PE_X | samtools view -Sb -F 4 - >$SAMSE_BAM"      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_sampe_sort"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA"
    HOLD="-hold_jid "$SAMPLE"_sampe"

    SAMPE_BAM=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_sampe.bam"
    SORTED_SAMPE=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_sorted_sampe.bam"
    
    CMD="${ENVIRONMENT[$ENV]}  samtools sort --threads $THREADS -o $SORTED_SAMPE $SAMPE_BAM "      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_samse_sort"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA"
    HOLD="-hold_jid "$SAMPLE"_samse"

    SAMSE_BAM=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_samse.bam"
    SORTED_SAMSE=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_sorted_samse.bam"
    
    CMD="${ENVIRONMENT[$ENV]} samtools sort --threads $THREADS -o $SORTED_SAMSE $SAMSE_BAM "      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_merge"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA"
    HOLD="-hold_jid "$SAMPLE"_samse_sort,"$SAMPLE"_sampe_sort"

    SORTED_SAMSE=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_sorted_samse.bam"
    SORTED_SAMPE=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_sorted_sampe.bam"
    MERGED_BAM=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE".bam"    

    CMD="${ENVIRONMENT[$ENV]} samtools merge $MERGED_BAM $SORTED_SAMSE $SORTED_SAMPE"      
    
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_add_rgs"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="SINGULARITY"
    HOLD="-hold_jid "$SAMPLE"_merge"
    
    MERGED_BAM=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE".bam"
    RG_BAM=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_rg.bam"        

    LANE=${RG_LANE[$SAMPLE]}
    CELL=${RG_CELL[$SAMPLE]}
    INDEX=${RG_INDEX[$SAMPLE]}

    CMD="${ENVIRONMENT[$ENV]} \
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

    ############
    JID=$SAMPLE"_mark_dups"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="SINGULARITY"
    HOLD="-hold_jid "$SAMPLE"_add_rgs"
    
    RG_BAM=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_rg.bam"        
    FINAL_BAM=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_processed.bam"        
    METRICS=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_processed.metrics"        


    CMD="${ENVIRONMENT[$ENV]} \
            gatk MarkDuplicates \
                --INPUT $RG_BAM \
                --OUTPUT $FINAL_BAM \
                --METRICS_FILE $METRICS \
                --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 900"
 
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_bam_index"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=1
    ENV="CONDA"
    HOLD="-hold_jid "$SAMPLE"_mark_dups"
    
    FINAL_BAM=$WORK_DIR"/results/01-processed_reads/outgroup_mapped_reads/"$SAMPLE"_processed.bam"        

    CMD="${ENVIRONMENT[$ENV]} samtools index $FINAL_BAM"
 
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"
   
done


