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

mkdir -p $SNP_DIR/r1

cd $SNP_DIR/r1

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
    OUT_VCF=$SNP_DIR/hc/$SAMPLE.vcf
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
#for INTERVAL in $(cat failed_1.intervals); do

    INTERVAL_DIR=$(echo $INTERVAL | sed 's/:/_/')

    ############ GDBIMPORT
    JID=$INTERVAL_DIR"_import_bsqr-1_fail-1"
    LOG=$LOGS_DIR/$JID".log"
    SCRIPT=$SCRIPTS_DIR/$JID".sh"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD=""

    SAMPLES_LIST=$SNP_DIR/samples.list
    OUT_DB=$SNP_DIR/db/$INTERVAL_DIR
    
    rm -r $OUT_DB

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
grep -L "Traversal complete" ../logs/*_import_bsqr-1.log >inc
#no failed jobs

#if failed
cat inc | cut -f3 -d/ | sed 's/_import_bsqr-1.log//' | sed 's/_/:/' >failed_1.intervals
#then re-run

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

    IN_DB=$SNP_DIR/db/$INTERVAL_DIR
    IN_REFERENCE=$HAE_GENOME
    OUT_VCF=$SNP_DIR/geno/$INTERVAL_DIR.vcf

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
ls $SNP_DIR/geno/*.vcf >$SNP_DIR/vcf.list

${ENVIRONMENTS['SINGULARITY']} \
    split -n l/600 \
        --numeric-suffixes=1 \
        --additional-suffix .list \
        $SNP_DIR/vcf.list \
        $SNP_DIR/merge-1_


for CHUNK in $(seq -w 1 600); do
    CMD="${ENVIRONMENTS["SINGULARITY"]} \
        gatk --java-options "-Xmx4g" MergeVcfs \
            --MAX_RECORDS_IN_RAM 500000 \
            -I merge-1_$CHUNK.list \
            -O merge-1_$CHUNK.vcf"

    echo $CMD | $QSUB -N merge-1_$CHUNK -o merge-1_$CHUNK.log -pe mpi 12
done

ls $SNP_DIR/merge-1_*.vcf >$SNP_DIR/vcf.list

${ENVIRONMENTS['SINGULARITY']} \
    split -n l/10 \
        --numeric-suffixes=1 \
        --additional-suffix .list \
        $SNP_DIR/vcf.list \
        $SNP_DIR/merge-2_

for CHUNK in $(seq -w 1 10); do
    CMD="${ENVIRONMENTS["SINGULARITY"]} \
        gatk --java-options "-Xmx4g" MergeVcfs \
            --MAX_RECORDS_IN_RAM 500000 \
            -I merge-2_$CHUNK.list \
            -O merge-2_$CHUNK.vcf"

    echo $CMD | $QSUB -N merge-2_$CHUNK -o merge-2_$CHUNK.log -pe mpi 12
done

#final merge (on titan)
ls $SNP_DIR/merge-2_*.vcf >$SNP_DIR/vcf.list

${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk --java-options "-Xmx24g" MergeVcfs \
        --MAX_RECORDS_IN_RAM 500000 \
        -I $SNP_DIR/vcf.list \
        -O $SNP_DIR/tmp.vcf

${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk --java-options "-Xmx24g" SortVcf \
        --MAX_RECORDS_IN_RAM 2000000 \
        -I $SNP_DIR/tmp.vcf \
        -O $SNP_DIR/cohort_raw_bqsr-1.vcf


#clean dir
rm merge-* tmp*

################################################################################

#filter snps
#select snps and filter
${ENVIRONMENTS["SINGULARITY"]} \
    gatk SelectVariants \
        -V $SNP_DIR/cohort_raw_bqsr-1.vcf \
        -select-type SNP \
        -O $SNP_DIR/filter/cohort_raw_bqsr-1_SNPs.vcf \
        -R $HAE_GENOME

${ENVIRONMENTS["SINGULARITY"]} \
    gatk SelectVariants \
        -V $SNP_DIR/cohort_raw_bqsr-1.vcf \
        -select-type INDEL \
        -O $SNP_DIR/filter/cohort_raw_bqsr-1_INDELs.vcf \
        -R $HAE_GENOME


############

${ENVIRONMENTS["SINGULARITY"]} \
    gatk VariantFiltration \
        -R $HAE_GENOME \
        -V $SNP_DIR/filter/cohort_raw_bqsr-1_SNPs.vcf \
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
        -O $SNP_DIR/filter/cohort_softFiltered_bqsr-1_SNPs.vcf

${ENVIRONMENTS["SINGULARITY"]} \
    gatk VariantFiltration \
            -R $HAE_GENOME \
            -V $SNP_DIR/filter/cohort_raw_bqsr-1_INDELs.vcf \
            --filter-name "QD_lt_2,indel" \
            --filter-expression "QD < 2.0" \
            --filter-name "FS_gt_200,indel" \
            --filter-expression "FS > 200.0" \
            --filter-name "ReadPosRankSum_lt_-20,indel" \
            --filter-expression "ReadPosRankSum < -20.0" \
            -O $SNP_DIR/filter/cohort_softFiltered_bqsr-1_INDELs.vcf

#and merge them back together.

ls $SNP_DIR/filter/cohort_softFiltered_bqsr-1_INDELs.vcf \
    $SNP_DIR/filter/cohort_softFiltered_bqsr-1_SNPs.vcf \
    >$SNP_DIR/filter/vcf.list

${ENVIRONMENTS["SINGULARITY"]} \
    gatk MergeVcfs \
        -I $SNP_DIR/filter/vcf.list \
        -O $SNP_DIR/filter/cohort_softFiltered_bqsr-1.vcf \
        -R $HAE_GENOME

vcftools \
    --remove-filtered-all \
    --vcf $SNP_DIR/filter/cohort_softFiltered_bqsr-1.vcf \
    --recode \
    --recode-INFO-all \
    --stdout \
    >$SNP_DIR/cohort_filtered_bqsr-1.vcf

#index
${ENVIRONMENTS["SINGULARITY"]} \
     gatk IndexFeatureFile \
        -F $SNP_DIR/cohort_filtered_bqsr-1.vcf


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
    IN_VCF=$SNP_DIR/cohort_filtered_bqsr-1.vcf
    OUT_TABLE=$SNP_DIR/recal/$SAMPLE"_bqsr-1_pre-cov.table"

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
    IN_VCF=$SNP_DIR/cohort_filtered_bqsr-1.vcf
    IN_TABLE=$SNP_DIR/recal/$SAMPLE"_bqsr-1_pre-cov.table"
    OUT_BAM=$SNP_DIR/recal/$SAMPLE"_bqsr-1.bam"

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

    IN_BAM=$SNP_DIR/recal/$SAMPLE"_bqsr-1.bam"
    IN_REFERENCE=$HAE_GENOME
    IN_VCF=$SNP_DIR/cohort_filtered_bqsr-1.vcf
    OUT_TABLE=$SNP_DIR/recal/$SAMPLE"_bqsr-1_post-cov.table"

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

    OUT_PDF=$SNP_DIR/recal/$SAMPLE"_bqsr-1.pdf"
    BEFORE_TABLE=$SNP_DIR/recal/$SAMPLE"_bqsr-1_pre-cov.table"
    AFTER_TABLE=$SNP_DIR/recal/$SAMPLE"_bqsr-1_post-cov.table"

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
ls "$SNP_DIR/recal/"*"_bqsr-1_pre-cov.table" >$SNP_DIR/recal/pre.list
ls "$SNP_DIR/recal/"*"_bqsr-1_post-cov.table" >$SNP_DIR/recal/post.list

${ENVIRONMENTS["SINGULARITY"]} \
    gatk GatherBQSRReports \
        --input $SNP_DIR/recal/pre.list \
        --output $SNP_DIR/recal/all_bqsr-1_pre-cov.table

${ENVIRONMENTS["SINGULARITY"]}
    gatk GatherBQSRReports \
        --input $SNP_DIR/recal/post.lis  \
        --output $SNP_DIR/recal/all_bqsr-1_post-cov.table
        
${ENVIRONMENTS["SINGULARITY"]}
gatk AnalyzeCovariates \
    -before $SNP_DIR/recal/all_bqsr-1_pre-cov.table \
    -after $SNP_DIR/recal/all_bqsr-1_post-cov.table \
    -plots $SNP_DIR/recal/all_bqsr-1.pdf


