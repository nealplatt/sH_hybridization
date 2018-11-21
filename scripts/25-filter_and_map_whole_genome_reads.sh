#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh

#START FILTERING/MAPPING PROCESS
cd $RESULTS_DIR

FILTER_DIR=$RESULTS_DIR/processed_reads/wgrs
MAP_DIR=$RESULTS_DIR/map_reads/wgrs
DATA_WGRS_DIR=$DATA_DIR/wgrs_data/AndersonT_08072018

mkdir -p $FILTER_DIR $MAP_DIR

#combine into single files per sample/read
for SAMPLE in "${WGRS_SAMPLES[@]}"; do

    PBS_FILE=$RESULTS_DIR/pbs_scripts/$SAMPLE"-filter_and_map.pbs.sh"
    rm $PBS_FILE

    echo '#!/bin/bash'                      >$PBS_FILE
    echo "#$ -V"                            >>$PBS_FILE
    echo "#$ -cwd"                          >>$PBS_FILE
    echo "#$ -S /bin/bash"                  >>$PBS_FILE
    echo "#$ -pe mpi 12"                    >>$PBS_FILE
    echo "#$ -j y"                          >>$PBS_FILE
    echo "#$ -q all.q"                      >>$PBS_FILE
    echo "#$ -N $SAMPLE.filter_map"         >>$PBS_FILE
    echo "#$ -o $SAMPLE.filter_map.log"     >>$PBS_FILE   
    echo -e "\n\n"                          >>$PBS_FILE 



    echo "#most of this work will be done in the conda environment" >>$PBS_FILE
    ENV="CONDA_SNPS"
    echo ${ENVIRONMENTS[$ENV]} >>$PBS_FILE
    echo "" >>$PBS_FILE 
    echo "" >>$PBS_FILE

    ############ ZCAT READS
    echo "#zcat data from seq files into a single file per read" >>$PBS_FILE
    for READ in R1 R2; do
        
        #output files
        OUT_FASTQ=$FILTER_DIR/$SAMPLE"_wgrs_raw_"$READ".fastq"

        CMD="#zcat $DATA_WGRS_DIR/"$SAMPLE"_S"*"_L00{1,2,3,4}_"$READ"_001.fastq.gz >"$OUT_FASTQ

        echo $CMD >>$PBS_FILE      
    done
    echo -e "\n\n" >>$PBS_FILE 

    ############ FILTER_READS
    echo "#filter reads with trimmomatic" >>$PBS_FILE   
    
    #input
    R1=$FILTER_DIR/$SAMPLE"_wgrs_raw_R1.fastq"
    R2=$FILTER_DIR/$SAMPLE"_wgrs_raw_R2.fastq"

    #output
    PE_R1=$FILTER_DIR/$SAMPLE"_wgrs_filtered_R1_paired.fastq.gz"
    PE_R2=$FILTER_DIR/$SAMPLE"_wgrs_filtered_R2_paired.fastq.gz"
    SE_R1=$FILTER_DIR/$SAMPLE"_wgrs_filtered_R1_unpaired.fastq.gz"
    SE_R2=$FILTER_DIR/$SAMPLE"_wgrs_filtered_R2_unpaired.fastq.gz"
    
    CMD="#trimmomatic PE -threads 12 -phred33 $R1 $R2 $PE_R1 $SE_R1 $PE_R2 $SE_R2 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15  MINLEN:36"
    
    echo $CMD >>$PBS_FILE
    echo -e "\n\n" >>$PBS_FILE      

    ############ ALN R1 and R2 READS
    echo "#use bwa to aln read to the genome" >>$PBS_FILE
    for READ in R1 R2; do

        #input
        FQ=$FILTER_DIR/$SAMPLE"_wgrs_filtered_"$READ"_paired.fastq.gz"
        GENOME=$HAE_GENOME

        #output
        SAI=$MAP_DIR/$SAMPLE"_"$READ".sai"

        CMD="#bwa aln -t 12 -f $SAI $GENOME $FQ"

        echo $CMD >>$PBS_FILE

    done
    echo -e "\n\n" >>$PBS_FILE      

    ############ ALN RX READS
    echo "#use bwa to aln read to the genome" >>$PBS_FILE

    #input
    GENOME=$HAE_GENOME
    SE_R1=$FILTER_DIR/$SAMPLE"_wgrs_filtered_R1_unpaired.fastq.gz"
    SE_R2=$FILTER_DIR/$SAMPLE"_wgrs_filtered_R2_unpaired.fastq.gz"
    SE=$FILTER_DIR/$SAMPLE"_wgrs_filtered_RX_unpaired.fastq.gz"
    #output
    SAI=$MAP_DIR/$SAMPLE"_RX.sai"

    CMD="#zcat $SE_R1 $SE_R2 >$SE; bwa aln -t 12 -f $SAI $GENOME $SE"

    echo $CMD >>$PBS_FILE
    echo -e "\n\n" >>$PBS_FILE      

    ############ SAMPE READS
    echo "#use bwa pair and map reads to genome" >>$PBS_FILE 
   
    #input
    GENOME=$HAE_GENOME
    SAI_1=$MAP_DIR/$SAMPLE"_R1.sai"
    SAI_2=$MAP_DIR/$SAMPLE"_R2.sai"
    PE_1=$FILTER_DIR/$SAMPLE"_wgrs_filtered_R1_paired.fastq.gz"
    PE_2=$FILTER_DIR/$SAMPLE"_wgrs_filtered_R2_paired.fastq.gz"

    #output
    SAMPE_BAM=$MAP_DIR/$SAMPLE"_wgrs_sampe.bam"

    CMD="#bwa sampe $GENOME $SAI_1 $SAI_2 $PE_1 $PE_2 | samtools view -Sb - >$SAMPE_BAM"

    echo $CMD >>$PBS_FILE
    echo -e "\n\n" >>$PBS_FILE      

    ############ SAMSE READS
    echo "#use bwa pair and map reads to genome" >>$PBS_FILE 
   
    #input
    GENOME=$HAE_GENOME
    SAI=$MAP_DIR/$SAMPLE"_RX.sai"
    SE=$FILTER_DIR/$SAMPLE"_wgrs_filtered_RX_unpaired.fastq.gz"

    #output
    SAMSE_BAM=$MAP_DIR/$SAMPLE"_wgrs_samse.bam"

    CMD="#bwa samse $GENOME $SAI $SE | samtools view -Sb - >$SAMSE_BAM"

    echo $CMD >>$PBS_FILE
    echo -e "\n\n" >>$PBS_FILE      


    ############ SORT SAMPE READS
    echo "#sort reads in sampe bam file" >>$PBS_FILE 
  
    #input
    SAMPE_BAM=$MAP_DIR/$SAMPLE"_wgrs_sampe.bam"

    #output
    SORTED_SAMPE=$MAP_DIR/$SAMPLE"_sorted_sampe.bam"

    CMD="#samtools sort --threads 10 -o $SORTED_SAMPE $SAMPE_BAM"

    echo $CMD >>$PBS_FILE
    echo -e "\n\n" >>$PBS_FILE      

    ############ SORT SAMSE READS
    echo "#sort reads in samse file" >>$PBS_FILE 
  
    #input
    SAMSE_BAM=$MAP_DIR/$SAMPLE"_wgrs_samse.bam"

    #output
    SORTED_SAMSE=$MAP_DIR/$SAMPLE"_sorted_samse.bam"

    CMD="#samtools sort --threads 10 -o $SORTED_SAMSE $SAMSE_BAM"

    echo $CMD >>$PBS_FILE
    echo -e "\n\n" >>$PBS_FILE      

    ############ MERGE SAMSE AND SAMSE
    echo "#merge sampe and samse bam files" >>$PBS_FILE 

    #in
    SORTED_SAMSE=$MAP_DIR/$SAMPLE"_sorted_samse.bam"
    SORTED_SAMPE=$MAP_DIR/$SAMPLE"_sorted_sampe.bam"
    
    #out
    MERGED_BAM=$MAP_DIR/$SAMPLE".bam"    

    CMD="#samtools merge $MERGED_BAM $SORTED_SAMSE $SORTED_SAMPE"      
    
    echo $CMD >>$PBS_FILE
    echo -e "\n\n" >>$PBS_FILE      

    ############ FLAGSTAT
    echo "#get basic mapping stats" >>$PBS_FILE 
   
    #input
    MERGED_BAM=$MAP_DIR/$SAMPLE".bam"    

    #output
    FLAG_STAT=$MAP_DIR/$SAMPLE".flagstat"

    CMD="#samtools flagstat $MERGED_BAM >$FLAG_STAT"

    echo $CMD >>$PBS_FILE
    echo -e "\n\n" >>$PBS_FILE      

    ############ FILTER UNMAPPED READS
    echo "#filter unmapped reads" >>$PBS_FILE 

    #input
    MERGED_BAM=$MAP_DIR/$SAMPLE".bam"    

    #output
    FILTERED_BAM=$MAP_DIR/$SAMPLE"_filtered.bam"

    CMD="#samtools view -@ 12 -b -F 4 $MERGED_BAM >$FILTERED_BAM"

    echo $CMD >>$PBS_FILE
    echo -e "\n\n" >>$PBS_FILE       


    ############ ADD READ GROUP INFO
    echo "#add readgroup info"          >>$PBS_FILE 
    echo "#! in singularity container"  >>$PBS_FILE 
    ENV="SINGULARITY"
    
    #input
    FILTERED_BAM=$MAP_DIR/$SAMPLE"_filtered.bam"

    #output    
    RG_BAM=$MAP_DIR/$SAMPLE"_rg.bam"        

    #PARAMETERS
    LANE="1"
    CELL="HWHGTBGX5"
    INDEX=${WGRS_RG_INDEX[$SAMPLE]}

    CMD="#${ENVIRONMENTS[$ENV]} gatk AddOrReplaceReadGroups \
                --INPUT=$FILTERED_BAM \
                --OUTPUT=$RG_BAM \
                --RGID=$CELL.$LANE \
                --RGLB=library1 \
                --RGPL=illumina \
                --RGPU=$CELL.$INDEX.$LANE  \
                --RGSM=$SAMPLE"

    echo $CMD >>$PBS_FILE
    echo -e "\n\n" >>$PBS_FILE      

    ############ MARK DUPLICATE READS
    echo "#identify duplicate reads"        >>$PBS_FILE 
    echo "#! in singularity container"      >>$PBS_FILE 
    ENV="SINGULARITY"

    #in
    RG_BAM=$MAP_DIR/$SAMPLE"_rg.bam"        
    
    #out
    FINAL_BAM=$MAP_DIR/$SAMPLE"_wgrs_processed.bam"        
    METRICS=$MAP_DIR/$SAMPLE"_wgrs_processed.metrics"        

    CMD="#${ENVIRONMENTS[$ENV]} \
            gatk MarkDuplicates \
                --INPUT $RG_BAM \
                --OUTPUT $FINAL_BAM \
                --METRICS_FILE $METRICS \
                --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 900"
 
    echo $CMD >>$PBS_FILE
    echo -e "\n\n" >>$PBS_FILE      

    ############ calculate genome coverage
    echo "#calculate genome coverage"       >>$PBS_FILE 

    #in
    FINAL_BAM=$MAP_DIR/$SAMPLE"_wgrs_processed.bam"        
    
    #out
    GENOME_COV=$MAP_DIR/$SAMPLE"_wgrs_processed.bam"        
    METRICS=$MAP_DIR/$SAMPLE"_wgrs_processed.metrics"        

    CMD="#${ENVIRONMENTS[$ENV]} \
            gatk MarkDuplicates \
                --INPUT $RG_BAM \
                --OUTPUT $FINAL_BAM \
                --METRICS_FILE $METRICS \
                --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 900"
 
    echo $CMD >>$PBS_FILE
    echo -e "\n\n" >>$PBS_FILE    

    ############ INDEX BAM FILE
    echo "#index bam file" >>$PBS_FILE 

    #in
    FINAL_BAM=$MAP_DIR/$SAMPLE"_wgrs_processed.bam"        

    #out
    GENOME_COV=$MAP_DIR/$SAMPLE"_genomecov".txt
    
    CMD="bedtools genomecov -ibam $FINAL_BAM -d >$COV_TXT"
 
    echo $CMD >>$PBS_FILE
    echo -e "\n\n" >>$PBS_FILE      

   
    ############ SUBMIT TO THE QUEUE
    qsub $PBS_FILE 
done




#the processed bam file *_wgrs_processed.bam will be used in downstream analyses
# including SNP calling and coverage calculations. 

################################################################################
################################################################################
################################################################################
#before running through GATK collect some basic stats that will be of use here
# and in other projects that use whole genome amplification (thinking spec of 
# elisha's work.

#first calculate raw genome coverage


#avg cov
for SAMPLE in "${WGRS_SAMPLES[@]}"; do
    echo -e $SAMPLE"\t"$(awk '{sum+=$3} END {print sum/NR}' $SAMPLE"_genomecov".txt) >>nuc_wgsr.cov &
done

#check mitochondria coverage
for SAMPLE in "${SAMPLES[@]}"; do
    echo -e $SAMPLE"\t"$(grep AMPZ01026399 $SAMPLE"_genome.cov" | awk '{sum+=$3} END {print sum/NR}') >>mito_wgsr.cov &
done 


#build coverage histogram
cut -f3 $SAMPLE"_schHae.genomeCov" | sort | uniq -c | awk 

#cov in non repetitive regions


#histogram of coverage (how much is to low to be used)

#how many of the non-mapping reads are derived from human

#compare to human
cd $DATA_DIR/genome
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
bwa index hg38.fa

cd $RESULTS_DIR/map_wgrs_reads

mkdir est_hg38_contam

cd est_hg38_contam

for SAMPLE in "${SAMPLES[@]}"; do
    for READ in R1 R2; do

        SAI=$SAMPLE"_"$READ"_hg38.sai"
        FQ=$RESULTS_DIR"/filter_wgrs_reads/"$SAMPLE"_wgrs_filtered_"$READ"_*paired.fastq.gz"
        GENOME=$DATA_DIR/genome/hg38.fa
        LOG=$SAMPLE"_"$READ"_bwa-aln-hg38.log"

        echo "source activate snp_calling; \
                bwa aln -t 12 -f $SAI $GENOME $FQ >$LOG 2>&1" \
            | $QSUB -pe mpi 12 -N $SAMPLE"_"$READ"_bwa-aln-hg38" -o $LOG

    done
done

for SAMPLE in "${SAMPLES[@]}"; do

    GENOME=$DATA_DIR/genome/hg38.fa
    SAI_1=$SAMPLE"_R1_hg38.sai"
    SAI_2=$SAMPLE"_R2_hg38.sai"
    PE_1=$RESULTS_DIR/filter_wgrs_reads/$SAMPLE"_wgrs_filtered_R1_paired.fastq.gz"
    PE_2=$RESULTS_DIR/filter_wgrs_reads/$SAMPLE"_wgrs_filtered_R2_paired.fastq.gz"
    SAMPE_BAM=$SAMPLE"_hg38.bam"
   
    #map PE
    echo "bwa sampe $GENOME $SAI_1 $SAI_2 $PE_1 $PE_2 \
        | samtools view -Sb - >$SAMPE_BAM" \
        |  $QSUB -pe mpi 12 -N $SAMPLE"_sampe_hg38" -o $SAMPLE"_sampe_hg38.log"

done

#get mapping stats
for SAMPLE in "${SAMPLES[@]}"; do

    SAMPE_BAM=$SAMPLE"_hg38.bam"

    samtools flagstat $SAMPE_BAM >$SAMPLE"_PE_hg38.flagstat" &

done

#sample, filt_reads, mapped_sH, mapped_human, sH_nuc_cov, sh_mito_cov
for SAMPLE in "${SAMPLES[@]}"; do
    #num reads
    READS=$(head -n1 $SAMPLE"_PE.flagstat" | cut -f1 -d" ")
    MAPPED=$(head -n5 $SAMPLE"_PE.flagstat" | tail -n 1 | cut -f1 -d" ")

    echo -e $SAMPLE"\t"$READS"\t"$MAPPED
done
    





   


    SAMPE_BAM=$SAMPLE"_wgrs_sampe.bam"

    echo "samtools flagstat $SAMPE_BAM >$SAMPLE.flagstat" \
        |  $QSUB -pe mpi 12 -N $SAMPLE"_flagstat" -o $SAMPLE"_flagstat.log" -hold_jid $SAMPLE"_sampe"






    #get mapping stats
    SAMPE_BAM=$SAMPLE"_wgrs_sampe.bam"

    echo "samtools flagstat $SAMPE_BAM >$SAMPLE.flagstat" \
        |  $QSUB -pe mpi 12 -N $SAMPLE"_flagstat" -o $SAMPLE"_flagstat.log" -hold_jid $SAMPLE"_sampe"

