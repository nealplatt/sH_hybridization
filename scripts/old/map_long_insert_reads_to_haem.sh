#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh

cd $FILTER_DIR

zcat ERR119622* | gzip >tmp_ERR119622.fastq.gz &
zcat ERR119623* | gzip >tmp_ERR119623.fastq.gz &
wait

mkdir long_insert_PE
mv ERR119622* long_insert_PE
mv ERR119623* long_insert_PE

mv tmp_ERR119622.fastq.gz ERR119622_filtered_unpaired_RX.fastq.gz
mv tmp_ERR119623.fastq.gz ERR119623_filtered_unpaired_RX.fastq.gz

cd $MAP_DIR

source activate snp_calling; bwa aln -t 12 -f ERR119622.sai $HAE_GENOME $FILTER_DIR/ERR119622_filtered_unpaired_RX.fastq.gz
bwa samse $HAE_GENOME ERR119622.sai_X $FILTER_DIR/ERR119622_filtered_unpaired_RX.fastq.gz | samtools view -Sb -F 4 - >ERR119622.bam      
samtools sort --threads 12 -o ERR119622_sorted.bam ERR119622.bam     
singularity exec /master/nplatt/schisto_hybridization/config/snpCalling_v0.0.8.img \
    gatk AddOrReplaceReadGroups \
        --INPUT=ERR119622_sorted.bam \
        --OUTPUT=ERR119622_rg.bam \
        --RGID=4.3 \
        --RGLB=library1 \
        --RGPL=illumina \
        --RGPU=4.ERR119622.3  \
        --RGSM=ERR119622
singularity exec /master/nplatt/schisto_hybridization/config/snpCalling_v0.0.8.img \
            gatk MarkDuplicates \
                --INPUT ERR119622_rg.bam \
                --OUTPUT ERR119622_processed.bam \
                --METRICS_FILE ERR119622_processed.metrics \
                --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 900
samtools index ERR119622_processed.bam



bwa aln -t 12 -f ERR119623.sai $HAE_GENOME $FILTER_DIR/ERR119623_filtered_unpaired_RX.fastq.gz
bwa samse $HAE_GENOME ERR119623.sai_X $FILTER_DIR/ERR119623_filtered_unpaired_RX.fastq.gz | samtools view -Sb -F 4 - >ERR119623.bam      
samtools sort --threads 12 -o ERR119623_sorted.bam ERR119623.bam     
singularity exec /master/nplatt/schisto_hybridization/config/snpCalling_v0.0.8.img \
    gatk AddOrReplaceReadGroups \
        --INPUT=ERR119623_sorted.bam \
        --OUTPUT=ERR119623_rg.bam \
        --RGID=5.2 \
        --RGLB=library1 \
        --RGPL=illumina \
        --RGPU=5.ERR119623.2  \
        --RGSM=ERR119623
singularity exec /master/nplatt/schisto_hybridization/config/snpCalling_v0.0.8.img \
            gatk MarkDuplicates \
                --INPUT ERR119623_rg.bam \
                --OUTPUT ERR119623_processed.bam \
                --METRICS_FILE ERR119623_processed.metrics \
                --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 900
samtools index ERR119623_processed.bam



