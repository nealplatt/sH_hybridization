#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh


#-------------------------------------------------------------------------------
#get a set of intervals that contain the probes:
cd $DATA_DIR

cp /master/fcheval/data/sh_exome/S0742423_Probes_w_mito.bed .

#these are overlapping probes so I need to combine them reads are less than 500bp
# so combining probes w/in 500bp of eo
bedtools merge -d 500 -i S0742423_Probes_w_mito.bed >merged_probes.bed

#from previous work I know that AMPZ01026399 is the mitochondrial contig...and
#DQ157222 was used in the probe design.  The real mito needs to be added and 
#ncbi contig removed.

grep AMPZ01026399 genome/schHae_v1.fa.fai
#AMPZ01026399.1  17526   355286760       80      81

sed '1 i\AMPZ01026399.1\t1\t17526' merged_probes.bed \
    | grep -v DQ157222 \
    >schHae_v1_probes.bed

awk '{print $1":"$2"-"$3}' schHae_v1_probes.bed >schHae_v1_probes.intervals
# this gives us 52,237 intervals that will be analyzed

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

mkdir -p $RECAL_DIR/r1

cd $RECAL_DIR/r1

mkdir hc db geno filter recal

#CALL SNPS AT REF_PANEL LOCI IN EACH OF THE SAMPLES (SH AND OUTGROUP)
for SAMPLE in "${SAMPLES[@]}"; do

    BAM=$MAP_DIR/$SAMPLE"_processed.bam"

    echo"" >$SCRIPTS_DIR/$SAMPLE"_hc.sh"

    ############ HC
    JID=$SAMPLE"_hc_bsqr-1"
    LOG=$LOGS_DIR/$JID".log"
    SCRIPT=$SCRIPTS_DIR/$JID".sh"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD=""

    IN_BAM=$MAP_DIR/$SAMPLE"_processed.bam"
    OUT_VCF=$RECAL_DIR/r1/hc/$SAMPLE.vcf
    IN_REFERENCE=$HAE_GENOME
    PROBE_LIST=$DATA_DIR/schHae_v1_probes.intervals

    
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
    JID=$INTERVAL_DIR"_import_bsqr-1"
    LOG=$LOGS_DIR/$JID".log"
    SCRIPT=$SCRIPTS_DIR/$JID".sh"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD=""

    SAMPLES_LIST=$RECAL_DIR/r1/samples.list
    OUT_DB=$RECAL_DIR/r1/db/$INTERVAL_DIR
    
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
grep -L "Traversal complete" ../../logs/*_import_bsqr-1.log >inc
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

    IN_DB=$RECAL_DIR/r1/db/$INTERVAL_DIR
    IN_REFERENCE=$HAE_GENOME
    OUT_VCF=$RECAL_DIR/r1/geno/$INTERVAL_DIR.vcf

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
ls $RECAL_DIR/r1/geno/*.vcf >$RECAL_DIR/r1/vcf.list

${ENVIRONMENTS['SINGULARITY']} \
    split -n l/600 \
        --numeric-suffixes=1 \
        --additional-suffix .list \
        $RECAL_DIR/r1/vcf.list \
        $RECAL_DIR/r1/merge-1_


for CHUNK in $(seq -w 1 600); do
    CMD="${ENVIRONMENTS["SINGULARITY"]} \
        gatk --java-options "-Xmx4g" MergeVcfs \
            --MAX_RECORDS_IN_RAM 500000 \
            -I merge-1_$CHUNK.list \
            -O merge-1_$CHUNK.vcf"

    echo $CMD | $QSUB -N merge-1_$CHUNK -o merge-1_$CHUNK.log -pe mpi 12
done

ls $RECAL_DIR/r1/merge-1_*.vcf >ls $RECAL_DIR/r1/vcf.list

${ENVIRONMENTS['SINGULARITY']} \
    split -n l/10 \
        --numeric-suffixes=1 \
        --additional-suffix .list \
        $RECAL_DIR/r1/vcf.list \
        $RECAL_DIR/r1/merge-2_

for CHUNK in $(seq -w 1 10); do
    CMD="${ENVIRONMENTS["SINGULARITY"]} \
        gatk --java-options "-Xmx4g" MergeVcfs \
            --MAX_RECORDS_IN_RAM 500000 \
            -I merge-2_$CHUNK.list \
            -O merge-2_$CHUNK.vcf"

    echo $CMD | $QSUB -N merge-2_$CHUNK -o merge-2_$CHUNK.log -pe mpi 12
done

ls $RECAL_DIR/r1/merge-1_*.vcf >ls $RECAL_DIR/r1/vcf.list

${ENVIRONMENTS["SINGULARITY"]} \
    gatk --java-options "-Xmx4g" MergeVcfs \
        --MAX_RECORDS_IN_RAM 500000 \
        -I $RECAL_DIR/r1/vcf.list \
        -O $RECAL_DIR/r1/tmp.vcf

${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk --java-options "-Xmx24g" SortVcf \
        --MAX_RECORDS_IN_RAM 500000 \
        -I $RECAL_DIR/r1/tmp.vcf \
        -O $RECAL_DIR/r1/cohort_raw_bqsr-1.vcf


#clean dir
rm merge-*

################################################################################

#filter snps
#select snps and filter
${ENVIRONMENTS["SINGULARITY"]} \
    gatk SelectVariants \
        -V $RECAL_DIR/r1/cohort_raw_bqsr-1.vcf \
        -select-type SNP \
        -O $RECAL_DIR/r1/filter/cohort_raw_bqsr-1_SNPs.vcf \
        -R $HAE_GENOME

${ENVIRONMENTS["SINGULARITY"]} \
    gatk SelectVariants \
        -V $RECAL_DIR/r1/cohort_raw_bqsr-1.vcf \
        -select-type INDEL \
        -O $RECAL_DIR/r1/filter/cohort_raw_bqsr-1_INDELs.vcf \
        -R $HAE_GENOME


############

${ENVIRONMENTS["SINGULARITY"]} \
    gatk VariantFiltration \
        -R $HAE_GENOME \
        -V $RECAL_DIR/r1/filter/cohort_raw_bqsr-1_SNPs.vcf \
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
        -O $RECAL_DIR/r1/filter/cohort_softFiltered_bqsr-1_SNPs.vcf

${ENVIRONMENTS["SINGULARITY"]} \
    gatk VariantFiltration \
            -R $HAE_GENOME \
            -V $RECAL_DIR/r1/filter/cohort_raw_bqsr-1_INDELs.vcf \
            --filter-name "QD_lt_2,indel" \
            --filter-expression "QD < 2.0" \
            --filter-name "FS_gt_200,indel" \
            --filter-expression "FS > 200.0" \
            --filter-name "ReadPosRankSum_lt_-20,indel" \
            --filter-expression "ReadPosRankSum < -20.0" \
            -O $RECAL_DIR/r1/filter/cohort_softFiltered_bqsr-1_INDELs.vcf

#and merge them back together.

ls $RECAL_DIR/r1/filter/cohort_softFiltered_bqsr-1_INDELs.vcf \
    $RECAL_DIR/r1/filter/cohort_softFiltered_bqsr-1_SNPs.vcf \
    >$RECAL_DIR/r1/filter/vcf.list

${ENVIRONMENTS["SINGULARITY"]} \
    gatk MergeVcfs \
        -I $RECAL_DIR/r1/filter/vcf.list \
        -O $RECAL_DIR/r1/filter/cohort_softFiltered_bqsr-1.vcf \
        -R $HAE_GENOME

vcftools \
    --remove-filtered-all \
    --vcf $RECAL_DIR/r1/filter/cohort_softFiltered_bqsr-1.vcf \
    --recode \
    --recode-INFO-all \
    --stdout \
    >$RECAL_DIR/r1/cohort_filtered_bqsr-1.vcf

#index
${ENVIRONMENTS["SINGULARITY"]} \
     gatk IndexFeatureFile \
        -F $RECAL_DIR/r1/cohort_filtered_bqsr-1.vcf


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
    IN_VCF=$RECAL_DIR/r1/cohort_filtered_bqsr-1.vcf
    OUT_TABLE=$RECAL_DIR/r1/recal/$SAMPLE"_bqsr-1_pre-cov.table"

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
    IN_VCF=$RECAL_DIR/r1/cohort_filtered_bqsr-1.vcf
    IN_TABLE=$RECAL_DIR/r1/recal/$SAMPLE"_bqsr-1_pre-cov.table"
    OUT_BAM=$RECAL_DIR/r1/recal/$SAMPLE"_bqsr-1.bam"

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

    IN_BAM=$RECAL_DIR/r1/recal/$SAMPLE"_bqsr-1.bam"
    IN_REFERENCE=$HAE_GENOME
    IN_VCF=$RECAL_DIR/r1/cohort_filtered_bqsr-1.vcf
    OUT_TABLE=$RECAL_DIR/r1/recal/$SAMPLE"_bqsr-1_post-cov.table"

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

    OUT_PDF=$RECAL_DIR/r1/recal/$SAMPLE"_bqsr-1.pdf"
    BEFORE_TABLE=$RECAL_DIR/r1/recal/$SAMPLE"_bqsr-1_pre-cov.table"
    AFTER_TABLE=$RECAL_DIR/r1/recal/$SAMPLE"_bqsr-1_post-cov.table"

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
ls "$RECAL_DIR/r1/recal/"*"_bqsr-1_pre-cov.table" >$RECAL_DIR/r1/recal/pre.list
ls "$RECAL_DIR/r1/recal/"*"_bqsr-1_post-cov.table" >$RECAL_DIR/r1/recal/post.list

${ENVIRONMENTS["SINGULARITY"]} \
    gatk GatherBQSRReports \
        --input $RECAL_DIR/r1/recal/pre.list \
        --output $RECAL_DIR/r1/recal/all_bqsr-1_pre-cov.table

${ENVIRONMENTS["SINGULARITY"]}
    gatk GatherBQSRReports \
        --input $RECAL_DIR/r1/recal/post.lis  \
        --output $RECAL_DIR/r1/recal/all_bqsr-1_post-cov.table
        
${ENVIRONMENTS["SINGULARITY"]}
gatk AnalyzeCovariates \
    -before $RECAL_DIR/r1/recal/all_bqsr-1_pre-cov.table \
    -after $RECAL_DIR/r1/recal/all_bqsr-1_post-cov.table \
    -plots $RECAL_DIR/r1/recal/all_bqsr-1.pdf


