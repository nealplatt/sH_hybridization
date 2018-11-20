#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh

cd $RESULTS_DIR
mkdir call_snps
cd call_snps

SAMPLES=(       "Sh.NE_Dai-002.1"       "Sh.NE_Dai-010.1" 
                "Sh.NE_Dai-013.3"       "Sh.NE_Dai-031.1"       "Sh.NE_Dai-033.1"
                "Sh.NE_Dai-044.1"       "Sh.NE_Dai-045.1"       "Sh.NE_Dai-051.1"   
                "Sh.NE_Dai-074.1"       "Sh.NE_Dai-146.1"       "Sh.NE_DaiCP-233.1"
                "Sh.NE_DaiCP-276.1"     "Sh.NE_DaiCP-277.2"     "Sh.NE_Doki-029.1" 
                "Sh.NE_Kar-001.1"       "Sh.NE_Kar-002.1"       "Sh.NE_Kar-076.1" 
                "Sh.NE_Kar-096.2"       "Sh.NE_Kar-241.1"       "Sh.NE_Kar-241.2" 
                "Sh.NE_Kar-281.1"       "Sh.NE_Kar-37.2"        "Sh.NE_Lata-007.3" 
                "Sh.NE_Lata-033.1"      "Sh.NE_Lata-078.1"      "Sh.NE_Lata-253.1" 
                "Sh.NE_Lata-275.2"      "Sh.NE_Lata-293.1"      "Sh.NE_Lata-294.1" 
                "Sh.NE_LibTB-009.2"     "Sh.NE_LibTB-010.1"     "Sh.NE_LibTB-022.1" 
                "Sh.NE_LibTB-028.1"     "Sh.NE_LibTB-031.1"     "Sh.NE_NG-011.1" 
                "Sh.NE_NG-06.2"         "Sh.NE_NG-089.1"        "Sh.NE_NG-236.1" 
                "Sh.NE_Seb-076.1"       "Sh.NE_Seb-078.2"       "Sh.NE_Seb-081.2" 
                "Sh.NE_Tiag-272.1"      "Sh.NE_YK-029.2"        "Sh.NE_YK-069.1" 
                "Sh.NE_YK-099.2"        "Sh.NE_YK-248.2"        "Sh.NE_Youri-069.2" 
                "Sh.NE_Youri-091.3"     "Sh.TZ_PEM0063.1"       "Sh.TZ_PEM0075.1" 
                "Sh.TZ_PEM0076.1"       "Sh.TZ_PEM0079.1"       "Sh.TZ_PEM0089.2" 
                "Sh.TZ_PEM0094.2"       "Sh.TZ_PEM0099.2"       "Sh.TZ_PEM0103.1" 
                "Sh.TZ_PEM0104.1"       "Sh.TZ_PEM0106.2"       "Sh.TZ_PEM0108.1" 
                "Sh.TZ_PEM0110.1"       "Sh.TZ_PEM0114.3"       "Sh.TZ_PEM0115.4" 
                "Sh.TZ_PEM0120.1"       "Sh.TZ_PEM0125.1"       "Sh.TZ_PEM0126.1" 
                "Sh.TZ_PEM0127.1"       "Sh.TZ_PEM0128.1"       "Sh.TZ_PEM0130.1" 
                "Sh.TZ_PEM0133.1"       "Sh.TZ_PEM0139.2"       "Sh.TZ_PEM0145.3" 
                "Sh.TZ_PEM0154.1"       "Sh.TZ_PEM0157.3"       "Sh.TZ_PEM0166.1" 
                "Sh.TZ_PEM0171.1"       "Sh.TZ_UNG0006.1"       "Sh.TZ_UNG0038.1"  
                "Sh.TZ_UNG0076.1"       "Sh.TZ_UNG0077.1"       "Sh.TZ_UNG0078.1" 
                "Sh.TZ_UNG0087.2"       "Sh.TZ_UNG0089.3"       "Sh.TZ_UNG0092.3" 
                "Sh.TZ_UNG0099.1"       "Sh.TZ_UNG0102.1"       "Sh.TZ_UNG0111.1" 
                "Sh.TZ_UNG0117.1"       "Sh.TZ_UNG0121.1"       "Sh.TZ_UNG0125.3" 
                "Sh.TZ_UNG0127.1"       "Sh.TZ_UNG0129.2"       "Sh.TZ_UNG0134.1" 
                "Sh.TZ_UNG0137.3"       "Sh.TZ_UNG0139.1"       "Sh.TZ_UNG0142.2" 
                "Sh.TZ_UNG0146.1"       "ERR119622"             "ERR103048" 
                "ERR119622"             "ERR310937"             "ERR119623" 
                "ERR084970"             "ERR037800"             "SRR433865" 
                "ERR103051"             "ERR119612"             "ERR119613" 
                "ERR310940"             "ERR539850"             "ERR539851" 
                "ERR539852"             "ERR539853"             "ERR539854" 
                "ERR539855"             "ERR539856"             "ERR539857"
                "Sh_Dai_044_1"          "Sh_DaiCP_276_1"        "Sh_Kar_001_1"
                "Sh_Kar_37_2"           "Sh_Lata_078_1"         "Sh_PEM_0103_1"
                "Sh_PEM_0104_1"         "Sh_PEM_0130_1"         "Sh_Tiag_272_1"
                "Sh_UNG_0038_1"         "Sh_UNG_0099_1"         "Sh_UNG_0121_1" )

mkdir hc db geno filter wgrs_db  wgrs_filter  wgrs_geno logs scripts 

############ HC
for SAMPLE in "${SAMPLES[@]}"; do

    JID="hc_"$SAMPLE
    LOG=$JID".log"
    THREADS=12
    ENV="SINGULARITY"     

    BAM=/master/nplatt/schisto_hybridization/results/map_reads/wgrs/$SAMPLE"_wgrs_processed.bam"
    OUT_VCF=./hc/$SAMPLE.vcf
    IN_REFERENCE=$HAE_GENOME
    
    CMD="${ENVIRONMENTS[$ENV]} \
        gatk HaplotypeCaller \
            -I $BAM \
            -O $OUT_VCF \
            -R $IN_REFERENCE \
            -ERC GVCF"

    echo $CMD | $QSUB -N hc_$SAMPLE -o hc_$SAMPLE.log -pe mpi 12

done

#WAIT UNTIL ALL JOBS ARE FINISHED.


#create list of samples
ls hc/*.vcf >samples.list


#run GDIMPORT for each contig
#get list of contigs in schHame genome
cut -f1 ../../data/genome/schHae_v1.fa.fai >contigs.list


#this needs to be run twice.  once with schisto exome and wgrs data
# and then once with only whole genome data

#WITH EXOMES
for INTERVAL in $(cat $RESULTS_DIR/call_snps/contigs.list); do

    ############ GDBIMPORT
    JID="import_"$INTERVAL
    LOG=$RESULTS_DIR/call_snps/logs/$JID".log"
    SCRIPT=$RESULTS_DIR/call_snps/scripts/$JID".sh"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD=""

    SAMPLES_LIST=$RESULTS_DIR/call_snps/samples.list
    OUT_DB=$RESULTS_DIR/call_snps/db/$INTERVAL
    
    if [ -d "$OUT_DB" ]; then
        rm -r $OUT_DB
        rm $LOG $SCRIPT
    fi 

    CMD="${ENVIRONMENTS[$ENV]} \
            gatk --java-options \"-Xmx4g -Xms4g\" GenomicsDBImport \
                -V $SAMPLES_LIST \
                --genomicsdb-workspace-path $OUT_DB \
                -L $INTERVAL \
                --reader-threads $THREADS \
                --batch-size 22"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"

    #limit num running jobs to 295
    NUM_JOBS=$(qstat | grep nplatt | wc -l)

    while [ $NUM_JOBS -gt 295 ]; do
        echo "."
        sleep 1
        NUM_JOBS=$(qstat | grep nplatt | wc -l)
    done

    echo $CMD | $JOB_QSUB
    echo $CMD >$SCRIPT
done

#check for completion
grep -L "Traversal complete" logs/import_*.log >inc
grep -i -e warn -e kill -e die -e error logs/import_*.log >>inc

#get a list of failed jobs to resubmit

############### ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

#then re-run
for INTERVAL in $(cat $RESULTS_DIR/call_snps/failed.list); do

    ############ GDBIMPORT
    JID="import_"$INTERVAL
    LOG=$RESULTS_DIR/call_snps/logs/$JID".log"
    SCRIPT=$RESULTS_DIR/call_snps/scripts/$JID".sh"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD=""

    SAMPLES_LIST=$RESULTS_DIR/call_snps/samples.list
    OUT_DB=$RESULTS_DIR/call_snps/db/$INTERVAL
    
    if [ -d "$OUT_DB" ]; then
        rm -r $OUT_DB
        rm $LOG $SCRIPT
    fi 

    CMD="${ENVIRONMENTS[$ENV]} \
            gatk --java-options \"-Xmx4g -Xms4g\" GenomicsDBImport \
                -V $SAMPLES_LIST \
                --genomicsdb-workspace-path $OUT_DB \
                -L $INTERVAL \
                --reader-threads $THREADS \
                --batch-size 22"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"

    echo $CMD | $JOB_QSUB
    echo $CMD >$SCRIPT
done


#clean up with tgs

#now import wgrs with outgroups (manually edited)
#cat wgrs_samples.list 
#hc/Sh_Dai_044_1.vcf
#hc/Sh_DaiCP_276_1.vcf
#hc/Sh_Kar_001_1.vcf
#hc/Sh_Kar_37_2.vcf
#hc/Sh_Lata_078_1.vcf
#hc/Sh_PEM_0103_1.vcf
#hc/Sh_PEM_0104_1.vcf
#hc/Sh_PEM_0130_1.vcf
#hc/Sh_Tiag_272_1.vcf
#hc/Sh_UNG_0038_1.vcf
#hc/Sh_UNG_0099_1.vcf
#hc/Sh_UNG_0121_1.vcf
#hc/ERR037800.vcf
#hc/ERR103048.vcf
#hc/ERR103051.vcf
#hc/ERR119612.vcf
#hc/ERR119613.vcf
#hc/ERR310937.vcf
#hc/ERR310940.vcf
#hc/ERR539855.vcf
#hc/ERR539857.vcf

for INTERVAL in $(cat $RESULTS_DIR/call_snps/contigs.list); do

    ############ GDBIMPORT
    JID="import_wgrs_"$INTERVAL
    LOG=$RESULTS_DIR/call_snps/logs/$JID".log"
    SCRIPT=$RESULTS_DIR/call_snps/scripts/$JID".sh"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD=""

    SAMPLES_LIST=$RESULTS_DIR/call_snps/wgrs_samples.list
    OUT_DB=$RESULTS_DIR/call_snps/wgrs_db/$INTERVAL
    
    if [ -d "$OUT_DB" ]; then
        rm -r $OUT_DB
        rm $LOG $SCRIPT
    fi 

    CMD="${ENVIRONMENTS[$ENV]} \
            gatk --java-options \"-Xmx4g -Xms4g\" GenomicsDBImport \
                -V $SAMPLES_LIST \
                --genomicsdb-workspace-path $OUT_DB \
                -L $INTERVAL \
                --reader-threads $THREADS \
                --batch-size 22"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"

    #limit num running jobs to 295
    NUM_JOBS=$(qstat | grep nplatt | wc -l)

    while [ $NUM_JOBS -gt 295 ]; do
        echo "."
        sleep 1
        NUM_JOBS=$(qstat | grep nplatt | wc -l)
    done

    echo $CMD | $JOB_QSUB
    echo $CMD >$SCRIPT
done

#check for completion
grep -L "Traversal complete" logs/import_*.log >inc
grep -i -e warn -e kill -e die -e error logs/import_*.log >>inc

#looks like 1535 and 1586 need to be "sleep"ed


#then genotype each interval for the exome + outroup data
for INTERVAL in $(cat $RESULTS_DIR/call_snps/contigs.list); do

    ############ genotype
    JID="geno_"$INTERVAL
    LOG=$RESULTS_DIR/call_snps/logs/$JID".log"
    SCRIPT=$RESULTS_DIR/call_snps/scripts/$JID".sh"
    THREADS=4
    ENV="SINGULARITY"     

    IN_DB=$RESULTS_DIR/call_snps/db/$INTERVAL
    IN_REFERENCE=$HAE_GENOME
    OUT_VCF=$SNP_DIR/geno/$INTERVAL.vcf

    CMD="${ENVIRONMENTS[$ENV]} \
            gatk GenotypeGVCFs \
                -R $IN_REFERENCE \
                -V gendb://$IN_DB \
                -new-qual \
                -O $OUT_VCF"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"

    #limit num running jobs
    NUM_JOBS=$(qstat | grep nplatt | wc -l)

    while [ $NUM_JOBS -gt 1000 ]; do
        echo "."
        sleep 1
        NUM_JOBS=$(qstat | grep nplatt | wc -l)
    done

    echo $CMD | $JOB_QSUB
    echo $CMD >$SCRIPT
done












for INTERVAL in $(cat $RESULTS_DIR/call_snps/contigs.list); do

    ############ genotype
    JID="geno_wgrs_"$INTERVAL
    LOG=$RESULTS_DIR/call_snps/logs/$JID".log"
    SCRIPT=$RESULTS_DIR/call_snps/scripts/$JID".sh"
    THREADS=4
    ENV="SINGULARITY"     

    IN_DB=$RESULTS_DIR/call_snps/wgrs_db/$INTERVAL
    IN_REFERENCE=$HAE_GENOME
    OUT_VCF=$SNP_DIR/geno/$INTERVAL.vcf

    CMD="${ENVIRONMENTS[$ENV]} \
            gatk GenotypeGVCFs \
                -R $IN_REFERENCE \
                -V gendb://$IN_DB \
                -new-qual \
                -O $OUT_VCF"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"

    #limit num running jobs
    NUM_JOBS=$(qstat | grep nplatt | wc -l)

    while [ $NUM_JOBS -gt 1000 ]; do
        echo "."
        sleep 1
        NUM_JOBS=$(qstat | grep nplatt | wc -l)
    done

    echo $CMD | $JOB_QSUB
    echo $CMD >$SCRIPT
done
























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


