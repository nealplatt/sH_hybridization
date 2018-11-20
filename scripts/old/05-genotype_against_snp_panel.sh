#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

#-------------------------------------------------------------------------------
#get a set of intervals that contain the probes:


mkdir $RESULTS_DIR/genotype_vs_panel

cd $RESULTS_DIR/genotype_vs_panel

mkdir hc db geno merge filter logs pbs

#CALL SNPS AT REF_PANEL LOCI IN EACH OF THE SAMPLES (SH AND OUTGROUP)
for SAMPLE in "${SAMPLES[@]}"; do

    BAM=$MAP_DIR/$SAMPLE"_processed.bam"

    echo"" >$SCRIPTS_DIR/$SAMPLE"_hc.sh"

    ############ HC
    JID=$SAMPLE"_hc_panel"
    LOG=logs/$JID".log"
    SCRIPT=pbs/$JID".sh"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD=""

    IN_BAM=$MAP_DIR/$SAMPLE"_processed.bam"
    OUT_VCF=hc/$SAMPLE.vcf
    IN_REFERENCE=$HAE_GENOME
    PROBE_LIST=$SNP_PANEL_DIR/snp_panel.list

    
    CMD="${ENVIRONMENTS[$ENV]} \
        gatk HaplotypeCaller \
            -I $BAM \
            -O $OUT_VCF \
            -R $IN_REFERENCE \
            -L $PROBE_LIST \
            -ERC GVCF"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >$SCRIPT
done

#WAIT UNTIL ALL JOBS ARE FINISHED.

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#create list of samples
ls hc/*.vcf >samples.list

#run GDIMPORT for each contig
for INTERVAL in $(cat $DATA_DIR/schHae_v1_probes.intervals); do

    INTERVAL_DIR=$(echo $INTERVAL | sed 's/:/_/')

    ############ GDBIMPORT
    JID=$INTERVAL_DIR"_import_panel"
    LOG=$LOGS_DIR/$JID".log"
    SCRIPT=$SCRIPTS_DIR/$JID".sh"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD=""

    SAMPLES_LIST=samples.list
    OUT_DB=db/$INTERVAL_DIR
    
    CMD="${ENVIRONMENTS[$ENV]} \
            gatk --java-options \"-Xmx4g -Xms4g\" GenomicsDBImport \
                -V $SAMPLES_LIST \
                --genomicsdb-workspace-path $OUT_DB \
                -L $INTERVAL \
                --reader-threads $THREADS \
                --batch-size 24"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    #limit num running jobs to 320
    NUM_JOBS=$(qstat | grep nplatt | wc -l)

    while [ $NUM_JOBS -gt 320 ]; do
        echo "."
        sleep 1
        NUM_JOBS=$(qstat | grep nplatt | wc -l)
    done

    echo $CMD | $JOB_QSUB
    echo $CMD >$SCRIPT
done

#check for completion
grep -L "Traversal complete" ../../logs/*_import_panel.log >inc
#no failed jobs

#then genotype each sample
for INTERVAL in $(cat $DATA_DIR/schHae_v1_probes.intervals); do

    INTERVAL_DIR=$(echo $INTERVAL | sed 's/:/_/')

    ############ genotype
    JID=$INTERVAL_DIR"_genotype-r1"
    LOG=$LOGS_DIR/$JID".log"
    SCRIPT=$SCRIPTS_DIR/$JID".sh"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD=""

    IN_DB=db/$INTERVAL_DIR
    IN_REFERENCE=$HAE_GENOME
    OUT_VCF=geno/$INTERVAL_DIR.vcf

    CMD="${ENVIRONMENTS[$ENV]} \
            gatk GenotypeGVCFs \
                -R $IN_REFERENCE \
                -V gendb://$IN_DB \
                -new-qual \
                -O $OUT_VCF"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    #limit num running jobs
    NUM_JOBS=$(qstat | grep nplatt | wc -l)

    while [ $NUM_JOBS -gt 500 ]; do
        echo "."
        sleep 1
        NUM_JOBS=$(qstat | grep nplatt | wc -l)
    done

    echo $CMD | $JOB_QSUB
    echo $CMD >$SCRIPT
done

#merge all contig vcf files into a single cohort file
ls geno/*.vcf >vcf.list

${ENVIRONMENTS['SINGULARITY']} \
    split -n l/600 \
        --numeric-suffixes=1 \
        --additional-suffix .list \
        vcf.list \
        merge-1_


for CHUNK in $(seq -w 1 600); do
    CMD="${ENVIRONMENTS["SINGULARITY"]} \
        gatk --java-options "-Xmx4g" MergeVcfs \
            --MAX_RECORDS_IN_RAM 500000 \
            -I merge-1_$CHUNK.list \
            -O merge-1_$CHUNK.vcf"

    echo $CMD | $QSUB -N merge-1_$CHUNK -o merge-1_$CHUNK.log -pe mpi 12
done

ls merge-1_*.vcf >ls vcf.list

${ENVIRONMENTS['SINGULARITY']} \
    split -n l/10 \
        --numeric-suffixes=1 \
        --additional-suffix .list \
        vcf.list \
        merge-2_

for CHUNK in $(seq -w 1 10); do
    CMD="${ENVIRONMENTS["SINGULARITY"]} \
        gatk --java-options "-Xmx4g" MergeVcfs \
            --MAX_RECORDS_IN_RAM 500000 \
            -I merge-2_$CHUNK.list \
            -O merge-2_$CHUNK.vcf"

    echo $CMD | $QSUB -N merge-2_$CHUNK -o merge-2_$CHUNK.log -pe mpi 12
done

ls merge-1_*.vcf >ls vcf.list

${ENVIRONMENTS["SINGULARITY"]} \
    gatk --java-options "-Xmx4g" MergeVcfs \
        --MAX_RECORDS_IN_RAM 500000 \
        -I vcf.list \
        -O tmp.vcf

${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk --java-options "-Xmx24g" SortVcf \
        --MAX_RECORDS_IN_RAM 500000 \
        -I tmp.vcf \
        -O cohort_raw_bqsr-1.vcf


#clean dir
rm merge-*

################################################################################

#filter snps
#select snps and filter
${ENVIRONMENTS["SINGULARITY"]} \
    gatk SelectVariants \
        -V cohort_raw_bqsr-1.vcf \
        -select-type SNP \
        -O filter/cohort_raw_bqsr-1_SNPs.vcf \
        -R $HAE_GENOME

${ENVIRONMENTS["SINGULARITY"]} \
    gatk SelectVariants \
        -V cohort_raw_bqsr-1.vcf \
        -select-type INDEL \
        -O filter/cohort_raw_bqsr-1_INDELs.vcf \
        -R $HAE_GENOME


############

${ENVIRONMENTS["SINGULARITY"]} \
    gatk VariantFiltration \
        -R $HAE_GENOME \
        -V filter/cohort_raw_bqsr-1_SNPs.vcf \
        --filter-name "QD_lt_2,snp" \
        --filter-expression "QD < 2.0" \
        --filter-name "FS_gt_60,snp" \
        --filter-expression "FS > 60.0" \
        --filter-name "MQ_lt_40,snp" \
        --filter-expression "MQ < 40.0" \
        --filter-name "MQRankSum_lt_-12.5,snp" \
        --filter-expression "MQRankSum < -12.5" \
        --filter-name "ReadPosRankSum_lt_-8,snp" \
        --filter-expression "ReadPosRankSum < -8.0" \
        -O filter/cohort_softFiltered_bqsr-1_SNPs.vcf

${ENVIRONMENTS["SINGULARITY"]} \
    gatk VariantFiltration \
            -R $HAE_GENOME \
            -V filter/cohort_raw_bqsr-1_INDELs.vcf \
            --filter-name "QD_lt_2,indel" \
            --filter-expression "QD < 2.0" \
            --filter-name "FS_gt_200,indel" \
            --filter-expression "FS > 200.0" \
            --filter-name "ReadPosRankSum_lt_-20,indel" \
            --filter-expression "ReadPosRankSum < -20.0" \
            -O filter/cohort_softFiltered_bqsr-1_INDELs.vcf

#and merge them back together.

ls filter/cohort_softFiltered_bqsr-1_INDELs.vcf \
    filter/cohort_softFiltered_bqsr-1_SNPs.vcf \
    >filter/vcf.list

${ENVIRONMENTS["SINGULARITY"]} \
    gatk MergeVcfs \
        -I filter/vcf.list \
        -O filter/cohort_softFiltered_bqsr-1.vcf \
        -R $HAE_GENOME

vcftools \
    --remove-filtered-all \
    --vcf filter/cohort_softFiltered_bqsr-1.vcf \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_filtered_bqsr-1.vcf

#index
${ENVIRONMENTS["SINGULARITY"]} \
     gatk IndexFeatureFile \
        -F cohort_filtered_bqsr-1.vcf


###################################
for SAMPLE in "${SAMPLES[@]}"; do

    BAM=$MAP_DIR/$SAMPLE"_processed.bam"

    ############ build table of covariates
    JID=$SAMPLE"_pre-covtable-r1"
    LOG=$LOGS_DIR/$JID".log"
    SCRIPT=$SCRIPTS_DIR/$JID".sh"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD=""

    IN_BAM=$BAM
    IN_REFERENCE=$HAE_GENOME
    IN_VCF=cohort_filtered_bqsr-1.vcf
    OUT_TABLE=recal/$SAMPLE"_bqsr-1_pre-cov.table"

    CMD="${ENVIRONMENTS["SINGULARITY"]} \
                gatk BaseRecalibrator \
                    --R $IN_REFERENCE \
                    -I $IN_BAM \
                    --known-sites $IN_VCF \
                    -O $OUT_TABLE"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"
    # with filtered cohort vcf...calculate covariates (per sample)
    
    echo $CMD | $JOB_QSUB
    echo $CMD >$SCRIPT

    ############ modify the quality scores in the bam files (per sample)
    JID=$SAMPLE"_modBam-r1"
    LOG=$LOGS_DIR/$JID".log"
    SCRIPT=$SCRIPTS_DIR/$JID".sh"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD="-hold_jid "$SAMPLE"_pre-covtable-r1"

    IN_BAM=$BAM
    IN_REFERENCE=$HAE_GENOME
    IN_VCF=cohort_filtered_bqsr-1.vcf
    IN_TABLE=recal/$SAMPLE"_bqsr-1_pre-cov.table"
    OUT_BAM=recal/$SAMPLE"_bqsr-1.bam"

    CMD="${ENVIRONMENTS["SINGULARITY"]} \
                gatk ApplyBQSR \
                    -R $IN_REFERENCE \
                    -I $IN_BAM \
                    --bqsr-recal-file $IN_TABLE \
                    -O $OUT_BAM"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"
    
    echo $CMD | $JOB_QSUB
    echo $CMD >$SCRIPT

    ############ #check for qual distribtuion after recal (per sample)
    JID=$SAMPLE"_post-covtable-r1"
    LOG=$LOGS_DIR/$JID".log"
    SCRIPT=$SCRIPTS_DIR/$JID".sh"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD="-hold_jid "$SAMPLE"_modBam-r1"

    IN_BAM=recal/$SAMPLE"_bqsr-1.bam"
    IN_REFERENCE=$HAE_GENOME
    IN_VCF=cohort_filtered_bqsr-1.vcf
    OUT_TABLE=recal/$SAMPLE"_bqsr-1_post-cov.table"

    CMD="${ENVIRONMENTS["SINGULARITY"]} \
                gatk BaseRecalibrator \
                    --R $IN_REFERENCE \
                    -I $IN_BAM \
                    --known-sites $IN_VCF \
                    -O $OUT_TABLE"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"
    
    echo $CMD | $JOB_QSUB
    echo $CMD >$SCRIPT

    ############# build plot to compare pre and post cal distributions (per sample)
    JID=$SAMPLE"_compare-r1"
    LOG=$LOGS_DIR/$JID".log"
    SCRIPT=$SCRIPTS_DIR/$JID".sh"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD="-hold_jid "$SAMPLE"_post-covtable-r1"

    OUT_PDF=recal/$SAMPLE"_bqsr-1.pdf"
    BEFORE_TABLE=recal/$SAMPLE"_bqsr-1_pre-cov.table"
    AFTER_TABLE=recal/$SAMPLE"_bqsr-1_post-cov.table"

    CMD="${ENVIRONMENTS["SINGULARITY"]} \
                gatk AnalyzeCovariates \
                -before $BEFORE_TABLE \
                -after $AFTER_TABLE \
                -plots $OUT_PDF"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"
    
    echo $CMD | $JOB_QSUB
    echo $CMD >$SCRIPT

done

# gather bqsr tables and examine for convergence across all samples
ls "recal/"*"_bqsr-1_pre-cov.table" >recal/pre.list
ls "recal/"*"_bqsr-1_post-cov.table" >recal/post.list

${ENVIRONMENTS["SINGULARITY"]} \
    gatk GatherBQSRReports \
        --input recal/pre.list \
        --output recal/all_bqsr-1_pre-cov.table

${ENVIRONMENTS["SINGULARITY"]}
    gatk GatherBQSRReports \
        --input recal/post.lis  \
        --output recal/all_bqsr-1_post-cov.table
        
${ENVIRONMENTS["SINGULARITY"]}
gatk AnalyzeCovariates \
    -before recal/all_bqsr-1_pre-cov.table \
    -after recal/all_bqsr-1_post-cov.table \
    -plots recal/all_bqsr-1.pdf


