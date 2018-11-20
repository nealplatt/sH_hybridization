declare -A ENVIRONMENTS
ENVIRONMENTS=(  ["SINGULARITY"]="singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img"
                ["CONDA_SNPS"]="source activate snp_calling;"
                ["CONDA_HYBRID_ID"]="source activate hybrid_id;"
                ["TITAN SINGULARITY"]="/opt/projects/singularity-2.4.2/bin/singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img"
             )

QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y"    

#set major directories
WORK_DIR="/master/nplatt/schisto_hybridization"
LOGS_DIR="logs"
RESULTS_DIR=$WORK_DIR"/results"
FILTER_DIR=$RESULTS_DIR"/01-processed_reads/filtered_reads"
MAP_DIR="./mt_panel"
DATA_DIR=$WORK_DIR"/data"
SEQ_DIR=$DATA_DIR"/sequence_data"
GENOME_DIR=$DATA_DIR"/genome"

#species id
mkdir logs

SAMPLES=(   
            "Sh.NE_YK-029.2"
            "Sh.NE_Dai-010.1"
            "Sh.NE_Kar-002.1"
            "Sh.NE_Doki-029.1"
            "ERR119623"
            "Sh.NE_Dai-051.1"
            "Sh.NE_Dai-146.1"
            "Sh.NE_DaiCP-276.1"
            "Sh.NE_Lata-275.2"
            "Sh.NE_Kar-241.1"
            "Sh.NE_Lata-078.1"
            "Sh.NE_Lata-294.1"
            "Sh.NE_LibTB-009.2"
            "Sh.NE_LibTB-010.1"
            "Sh.TZ_PEM0063.1"
            "Sh.TZ_PEM0094.2"
            "Sh.TZ_PEM0099.2"
            "Sh.TZ_PEM0103.1"
            "Sh.TZ_PEM0104.1"
            "Sh.TZ_PEM0106.2"
            "Sh.TZ_PEM0110.1"
            "Sh.TZ_PEM0114.3"
            "Sh.TZ_PEM0115.4"
            "Sh.TZ_PEM0120.1"
            "Sh.TZ_PEM0125.1"
            "Sh.TZ_PEM0126.1"
            "Sh.TZ_PEM0128.1"
            "Sh.TZ_PEM0130.1"
            "Sh.TZ_PEM0133.1"
            "Sh.TZ_PEM0145.3"
            "Sh.TZ_PEM0154.1"
            "Sh.TZ_PEM0157.3"
            "Sh.TZ_PEM0166.1"
            "Sh.TZ_PEM0171.1"
            "Sh.TZ_UNG0038.1"
            "Sh.TZ_UNG0076.1"
            "Sh.TZ_UNG0077.1"
            "Sh.TZ_UNG0092.3"
            "Sh.TZ_UNG0117.1"
            "Sh.TZ_UNG0137.3"
            "Sh.TZ_UNG0146.1"
            "Sh.NE_Dai-033.1"
            "Sh.TZ_UNG0129.2"
            "ERR103048"
            "ERR119622"
            "Sh.NE_YK-099.2"
            "Sh.NE_YK-069.1"
            "Sh.NE_YK-248.2"
            "Sh.NE_DaiCP-277.2"
            "Sh.NE_Dai-074.1"
            "Sh.NE_LibTB-031.1"
            "Sh.NE_LibTB-028.1"
            "Sh.NE_LibTB-022.1"
            "Sh.NE_NG-089.1"
            "Sh.NE_Youri-069.2"
            "Sh.NE_Kar-076.1"
            "Sh.NE_Kar-281.1"
            "Sh.NE_Kar-001.1"
            "Sh.NE_Seb-078.2"
            "Sh.NE_Seb-076.1"
            "Sh.NE_Seb-081.2"
            "Sh.TZ_UNG0099.1"
            "Sh.TZ_UNG0102.1"
            "Sh.TZ_UNG0111.1"
            "Sh.TZ_UNG0121.1"
            "Sh.TZ_UNG0125.3"
            "Sh.TZ_UNG0127.1"
            "Sh.TZ_UNG0134.1"
            "Sh.TZ_UNG0139.1"
            "ERR037800"
            "ERR084970"
            "SRR433865"
            "Sh.NE_Kar-37.2"
            "Sh.NE_Tiag-272.1"
            "Sh.NE_Lata-007.3"
            "Sh.NE_Kar-096.2"
            "Sh.NE_Lata-033.1"
            "Sh.NE_Kar-241.2"
            "Sh.NE_Dai-002.1"
            "Sh.NE_NG-06.2"
            "Sh.NE_Dai-013.3"
            "Sh.NE_Dai-031.1"
            "Sh.NE_Lata-253.1"
            "Sh.NE_Dai-044.1"
            "Sh.NE_NG-011.1"
            "Sh.NE_Dai-045.1"
            "Sh.NE_NG-236.1"
            "Sh.NE_Lata-293.1"
            "Sh.NE_Youri-091.3"
            "Sh.NE_DaiCP-233.1"
            "Sh.TZ_UNG0142.2"
            "Sh.TZ_PEM0076.1"
            "Sh.TZ_PEM0079.1"
            "Sh.TZ_PEM0089.2"
            "Sh.TZ_PEM0108.1"
            "Sh.TZ_PEM0127.1"
            "Sh.TZ_PEM0139.2"
            "Sh.TZ_UNG0006.1"
            "Sh.TZ_UNG0078.1"
            "Sh.TZ_UNG0087.2"
            "Sh.TZ_UNG0089.3"
)
#bowtie to reference panel
for SAMPLE in "${SAMPLES[@]}"; do

    HOLD=""
    echo"" >$SAMPLE"_process.sh"

    ############ MAP READS
    for READ in R1 R2; do
        JID=$SAMPLE"_"$READ"_map"
        LOG=$LOGS_DIR/$JID".log"
        THREADS=12 
        ENV="CONDA_SNPS"
    
        SAI=$MAP_DIR/$SAMPLE"_"$READ".sai"
        FQ=$FILTER_DIR/$SAMPLE"_filtered_paired_"$READ".fastq.gz"
        GENOME=schistosoma_mt_references.fasta

        JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"

        CMD="${ENVIRONMENTS[$ENV]} bwa aln -t 12 -f $SAI $GENOME $FQ"
            
        echo $CMD | $JOB_QSUB
        echo $CMD >>$SAMPLE"_process.sh"
    done

    ############ MAP READS
    JID=$SAMPLE"_RX_map"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12 
    ENV="CONDA_SNPS"
    
    SAI=$MAP_DIR/$SAMPLE"_RX.sai"
    FQ=$FILTER_DIR/$SAMPLE"_filtered_unpaired_RX.fastq.gz"
    GENOME=schistosoma_mt_references.fasta

    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"

    CMD="${ENVIRONMENTS[$ENV]} bwa aln -t 12 -f $SAI $GENOME $FQ"
            
    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_sampe"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"     
    HOLD="-hold_jid "$SAMPLE"_R1_map,"$SAMPLE"_R2_map"

    GENOME=schistosoma_mt_references.fasta
    SAI_1=$MAP_DIR/$SAMPLE"_R1.sai"
    SAI_2=$MAP_DIR/$SAMPLE"_R2.sai"
    PE_1=$FILTER_DIR/$SAMPLE"_filtered_paired_R1.fastq.gz"
    PE_2=$FILTER_DIR/$SAMPLE"_filtered_paired_R2.fastq.gz"
    SAMPE_BAM=$MAP_DIR/$SAMPLE"_sampe.bam"
    
    CMD="${ENVIRONMENTS[$ENV]} bwa sampe $GENOME $SAI_1 $SAI_2 $PE_1 $PE_2 | samtools view -Sb -F 4 - >$SAMPE_BAM"      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_samse"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_RX_map"

    GENOME=schistosoma_mt_references.fasta
    SAI_X=$MAP_DIR/$SAMPLE"_RX.sai"
    PE_X=$FILTER_DIR/$SAMPLE"_filtered_unpaired_RX.fastq.gz"
    SAMSE_BAM=$MAP_DIR/$SAMPLE"_samse.bam"
    
    CMD="${ENVIRONMENTS[$ENV]} bwa samse $GENOME $SAI_X $PE_X | samtools view -Sb -F 4 - >$SAMSE_BAM"      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_sampe_sort"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_sampe"

    SAMPE_BAM=$MAP_DIR/$SAMPLE"_sampe.bam"
    SORTED_SAMPE=$MAP_DIR/$SAMPLE"_sorted_sampe.bam"
    
    CMD="${ENVIRONMENTS[$ENV]}  samtools sort --threads $THREADS -o $SORTED_SAMPE $SAMPE_BAM "      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
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

    ############
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
done


rm *.sai
rm *_sam?e.bam

bwa index curassoni_mt_reference.fasta 
bwa index haematobium_mt_reference.fasta 


MAP_DIR=./haem_mito
LOGS_DIR=./haem_mito/logs
mkdir $MAP_DIR $LOGS_DIR

#bowtie to haem mito 
for SAMPLE in "${SAMPLES[@]}"; do

    HOLD=""
    echo"" >$SAMPLE"_process.sh"

    ############ MAP READS
    for READ in R1 R2; do
        JID=$SAMPLE"_"$READ"_haem_map"
        LOG=$LOGS_DIR/$JID".log"
        THREADS=12 
        ENV="CONDA_SNPS"
    
        SAI=$MAP_DIR/$SAMPLE"_"$READ".sai"
        FQ=$FILTER_DIR/$SAMPLE"_filtered_paired_"$READ".fastq.gz"
        GENOME=haematobium_mt_reference.fasta

        JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"

        CMD="${ENVIRONMENTS[$ENV]} bwa aln -t 12 -f $SAI $GENOME $FQ"
            
        echo $CMD | $JOB_QSUB
        echo $CMD >>$SAMPLE"_process.sh"
    done

    ############ MAP READS
    JID=$SAMPLE"_RX_haem_map"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12 
    ENV="CONDA_SNPS"
    
    SAI=$MAP_DIR/$SAMPLE"_RX.sai"
    FQ=$FILTER_DIR/$SAMPLE"_filtered_unpaired_RX.fastq.gz"
    GENOME=haematobium_mt_reference.fasta

    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"

    CMD="${ENVIRONMENTS[$ENV]} bwa aln -t 12 -f $SAI $GENOME $FQ"
            
    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_sampe_haem_"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"     
    HOLD="-hold_jid "$SAMPLE"_R1_haem_map,"$SAMPLE"_R2_haem_map"

    GENOME=haematobium_mt_reference.fasta
    SAI_1=$MAP_DIR/$SAMPLE"_R1.sai"
    SAI_2=$MAP_DIR/$SAMPLE"_R2.sai"
    PE_1=$FILTER_DIR/$SAMPLE"_filtered_paired_R1.fastq.gz"
    PE_2=$FILTER_DIR/$SAMPLE"_filtered_paired_R2.fastq.gz"
    SAMPE_BAM=$MAP_DIR/$SAMPLE"_sampe.bam"
    
    CMD="${ENVIRONMENTS[$ENV]} bwa sampe $GENOME $SAI_1 $SAI_2 $PE_1 $PE_2 | samtools view -Sb -F 4 - >$SAMPE_BAM"      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_samse_haem_"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_RX_map_haem_"

    GENOME=haematobium_mt_reference.fasta
    SAI_X=$MAP_DIR/$SAMPLE"_RX.sai"
    PE_X=$FILTER_DIR/$SAMPLE"_filtered_unpaired_RX.fastq.gz"
    SAMSE_BAM=$MAP_DIR/$SAMPLE"_samse.bam"
    
    CMD="${ENVIRONMENTS[$ENV]} bwa samse $GENOME $SAI_X $PE_X | samtools view -Sb -F 4 - >$SAMSE_BAM"      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_sampe_sort_haem_"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_sampe_haem_"

    SAMPE_BAM=$MAP_DIR/$SAMPLE"_sampe.bam"
    SORTED_SAMPE=$MAP_DIR/$SAMPLE"_sorted_sampe.bam"
    
    CMD="${ENVIRONMENTS[$ENV]}  samtools sort --threads $THREADS -o $SORTED_SAMPE $SAMPE_BAM "      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_samse_sort_haem_"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_samse_haem_"

    SAMSE_BAM=$MAP_DIR/$SAMPLE"_samse.bam"
    SORTED_SAMSE=$MAP_DIR/$SAMPLE"_sorted_samse.bam"
    
    CMD="${ENVIRONMENTS[$ENV]} samtools sort --threads $THREADS -o $SORTED_SAMSE $SAMSE_BAM "      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_merge_haem_"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_samse_sort_haem_,"$SAMPLE"_sampe_sort_haem_"

    SORTED_SAMSE=$MAP_DIR/$SAMPLE"_sorted_samse.bam"
    SORTED_SAMPE=$MAP_DIR/$SAMPLE"_sorted_sampe.bam"
    MERGED_BAM=$MAP_DIR/$SAMPLE".bam"    

    CMD="${ENVIRONMENTS[$ENV]} samtools merge $MERGED_BAM $SORTED_SAMSE $SORTED_SAMPE"      
    
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"
done

MAP_DIR=./curs_mito
LOGS_DIR=./curs_mito/logs
mkdir $MAP_DIR $LOGS_DIR

#bowtie to curs mito 
for SAMPLE in "${SAMPLES[@]}"; do

    HOLD=""
    echo"" >$SAMPLE"_process.sh"

    ############ MAP READS
    for READ in R1 R2; do
        JID=$SAMPLE"_"$READ"_curs_map"
        LOG=$LOGS_DIR/$JID".log"
        THREADS=12 
        ENV="CONDA_SNPS"
    
        SAI=$MAP_DIR/$SAMPLE"_"$READ".sai"
        FQ=$FILTER_DIR/$SAMPLE"_filtered_paired_"$READ".fastq.gz"
        GENOME=curassoni_mt_reference.fasta

        JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"

        CMD="${ENVIRONMENTS[$ENV]} bwa aln -t 12 -f $SAI $GENOME $FQ"
            
        echo $CMD | $JOB_QSUB
        echo $CMD >>$SAMPLE"_process.sh"
    done

    ############ MAP READS
    JID=$SAMPLE"_RX_curs_map"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12 
    ENV="CONDA_SNPS"
    
    SAI=$MAP_DIR/$SAMPLE"_RX.sai"
    FQ=$FILTER_DIR/$SAMPLE"_filtered_unpaired_RX.fastq.gz"
    GENOME=curassoni_mt_reference.fasta

    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"

    CMD="${ENVIRONMENTS[$ENV]} bwa aln -t 12 -f $SAI $GENOME $FQ"
            
    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_sampe_curs_"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"     
    HOLD="-hold_jid "$SAMPLE"_R1_curs_map,"$SAMPLE"_R2_curs_map"

    GENOME=curassoni_mt_reference.fasta
    SAI_1=$MAP_DIR/$SAMPLE"_R1.sai"
    SAI_2=$MAP_DIR/$SAMPLE"_R2.sai"
    PE_1=$FILTER_DIR/$SAMPLE"_filtered_paired_R1.fastq.gz"
    PE_2=$FILTER_DIR/$SAMPLE"_filtered_paired_R2.fastq.gz"
    SAMPE_BAM=$MAP_DIR/$SAMPLE"_sampe.bam"
    
    CMD="${ENVIRONMENTS[$ENV]} bwa sampe $GENOME $SAI_1 $SAI_2 $PE_1 $PE_2 | samtools view -Sb -F 4 - >$SAMPE_BAM"      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_samse_curs_"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_RX_map_curs_"

    GENOME=curassoni_mt_reference.fasta
    SAI_X=$MAP_DIR/$SAMPLE"_RX.sai"
    PE_X=$FILTER_DIR/$SAMPLE"_filtered_unpaired_RX.fastq.gz"
    SAMSE_BAM=$MAP_DIR/$SAMPLE"_samse.bam"
    
    CMD="${ENVIRONMENTS[$ENV]} bwa samse $GENOME $SAI_X $PE_X | samtools view -Sb -F 4 - >$SAMSE_BAM"      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_sampe_sort_curs_"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_sampe_curs_"

    SAMPE_BAM=$MAP_DIR/$SAMPLE"_sampe.bam"
    SORTED_SAMPE=$MAP_DIR/$SAMPLE"_sorted_sampe.bam"
    
    CMD="${ENVIRONMENTS[$ENV]}  samtools sort --threads $THREADS -o $SORTED_SAMPE $SAMPE_BAM "      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_samse_sort_curs_"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_samse_curs_"

    SAMSE_BAM=$MAP_DIR/$SAMPLE"_samse.bam"
    SORTED_SAMSE=$MAP_DIR/$SAMPLE"_sorted_samse.bam"
    
    CMD="${ENVIRONMENTS[$ENV]} samtools sort --threads $THREADS -o $SORTED_SAMSE $SAMSE_BAM "      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"

    ############
    JID=$SAMPLE"_merge_curs_"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="CONDA_SNPS"
    HOLD="-hold_jid "$SAMPLE"_samse_sort_curs_,"$SAMPLE"_sampe_sort_curs_"

    SORTED_SAMSE=$MAP_DIR/$SAMPLE"_sorted_samse.bam"
    SORTED_SAMPE=$MAP_DIR/$SAMPLE"_sorted_sampe.bam"
    MERGED_BAM=$MAP_DIR/$SAMPLE".bam"    

    CMD="${ENVIRONMENTS[$ENV]} samtools merge $MERGED_BAM $SORTED_SAMSE $SORTED_SAMPE"      
    
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SAMPLE"_process.sh"
done

echo -e "SAMPLE\tnum_raw\tnum_filtered\tnum_mapped(SHg)" >counts.table
for SAMPLE in "${SAMPLES[@]}"; do
    RAW_R1="/master/nplatt/schisto_hybridization/data/sequence_data/"$SAMPLE"_R1.fastq.gz"
    RAW_R2="/master/nplatt/schisto_hybridization/data/sequence_data/"$SAMPLE"_R2.fastq.gz"
    FILTER_R1=$FILTER_DIR/$SAMPLE"_filtered_paired_R1.fastq.gz"
    FILTER_R2=$FILTER_DIR/$SAMPLE"_filtered_paired_R2.fastq.gz"
    FILTER_RX=$FILTER_DIR/$SAMPLE"_filtered_unpaired_RX.fastq.gz"

    RAW_R1_COUNTS=$(($(zcat $RAW_R1 | wc -l) / 4)) 
    RAW_R2_COUNTS=$(($(zcat $RAW_R2 | wc -l) / 4)) 
    FILTER_R1_COUNTS=$(($(zcat $FILTER_R1 | wc -l) / 4)) 
    FILTER_R2_COUNTS=$(($(zcat $FILTER_R2 | wc -l) / 4)) 
    FILTER_RX_COUNTS=$(($(zcat $FILTER_RX | wc -l) / 4)) 

    MAPPED_READS=$(samtools view /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/$SAMPLE"_processed.bam" | wc -l)
    
    TOTAL_RAW=$((RAW_R1_COUNTS + RAW_R2_COUNTS))
    TOTAL_FILTER=$((FILTER_R1_COUNTS + FILTER_R2_COUNTS + FILTER_RX_COUNTS))
    TOTAL_MAP=$MAPPED_READS


    echo -e "$SAMPLE $TOTAL_RAW $TOTAL_FILTER $TOTAL_MAP" >>counts.table
done

mkdir int_files
mkdir int_files
mv Sh*sai int_files/
mv Sh*sam*bam int_files/


#histogram of quality scores (per sample)

#histogram of length (per sample)
