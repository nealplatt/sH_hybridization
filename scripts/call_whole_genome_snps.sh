#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR
mkdir call_whole_genome_snps
cd call_whole_genome_snps

mkdir hc db geno filter logs scripts 

SAMPLES=(   "Sh_Dai_044_1"          "Sh_DaiCP_276_1"        "Sh_Kar_001_1"
            "Sh_Kar_37_2"           "Sh_Lata_078_1"         "Sh_PEM_0103_1"
            "Sh_PEM_0104_1"         "Sh_PEM_0130_1"         "Sh_Tiag_272_1"
            "Sh_UNG_0038_1"         "Sh_UNG_0099_1"         "Sh_UNG_0121_1" )


############ HC
for SAMPLE in "${SAMPLES[@]}"; do

    #properties
    JID="hc_"$SAMPLE
    LOG=$JID".log"
    THREADS=12

    #input
    BAM=/master/nplatt/schisto_hybridization/results/map_reads/wgrs/$SAMPLE"_wgrs_processed.bam"
    IN_REFERENCE=$HAE_GENOME

    #output
    OUT_VCF=./hc/$SAMPLE-whole_genome.vcf

    #command
    /opt/projects/singularity-2.4.2/bin/singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
        gatk --java-options "-Xmx4g -Xms4g" HaplotypeCaller \
            -I $BAM \
            -O $OUT_VCF \
            -R $IN_REFERENCE \
            -ERC GVCF >$SAMPLE-hc.log 2>&1 &

#    CMD="source activate snp_calling; \
#            singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
#                gatk HaplotypeCaller \
#                    -I $BAM \
#                    -O $OUT_VCF \
#                    -R $IN_REFERENCE \
#                    -ERC GVCF"

    
# #submit
#    echo $CMD | $QSUB -N $JID -o $JID.log -pe mpi $THREADS

done

#recall snps on bovis whole genome data
#input
BAM=/master/nplatt/schisto_hybridization/results/map_reads/exome/ERR103048_processed.bam
IN_REFERENCE=$HAE_GENOME

#output
OUT_VCF=$(pwd)/hc/ERR103048-whole_genome.vcf

#command
/opt/projects/singularity-2.4.2/bin/singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    gatk --java-options "-Xmx4g -Xms4g" HaplotypeCaller \
        -I $BAM \
        -O $OUT_VCF \
        -R $IN_REFERENCE \
        -ERC GVCF \
        >ERR103048-hc.log 2>&1 &


#move log files
mv *.log logs

#-------------------------------------------------------------------------------

#create list of samples
ls hc/*.vcf >samples.list

#run GDIMPORT for each contig
#get list of contigs in schHame genome
cut -f1 ../../data/genome/schHae_v1.fa.fai >contigs.list

################################################################################

#import data
rm -r db/*
for INTERVAL in $(cat contigs.list); do

    ############ GDBIMPORT
    JID="import_"$INTERVAL
    LOG=$(pwd)/logs/$JID.log


    SAMPLES_LIST=$(pwd)/samples.list
    OUT_DB=$(pwd)/db/$INTERVAL
    

    CMD="singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
            gatk --java-options \"-Xmx4g -Xms4g\" GenomicsDBImport \
                -V $SAMPLES_LIST \
                --genomicsdb-workspace-path $OUT_DB \
                -L $INTERVAL \
                --reader-threads 4 \
                --batch-size 13"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi 12"

    #limit num running jobs to 600
    NUM_JOBS=$(qstat | grep nplatt | wc -l)

    while [ $NUM_JOBS -gt 600 ]; do
        echo "."
        sleep 1
        NUM_JOBS=$(qstat | grep nplatt | wc -l)
    done

    echo $CMD | $JOB_QSUB

done

#check for completion
rm failed.list

for INTERVAL in $(cat contigs.list); do
    if [ ! -d ./db/$INTERVAL ]; then
        echo $INTERVAL >>failed.list
    fi
done

#check for completion, get a list of failed jobs to resubmit
grep -L "Traversal complete" logs/import_*.log >>failed.list
grep -i -e warn -e kill -e die logs/import_*.log >>failed.list
sed -i 's/logs\/import_\(.*\).log/\1/gi' failed.list
sort -u failed.list >tmp; mv tmp failed.list

for INTERVAL in $(cat ./failed.list); do

    ############ GDBIMPORT
    JID="import_"$INTERVAL
    OUT_DB=$(pwd)/db/$INTERVAL
    
    #delete db from previous run
    if [ -d "$OUT_DB" ]; then
        rm -r $OUT_DB
    fi

    #delete log from previous run
    if [ -f logs/$JID.log ]; then
        rm logs/$JID.log
    fi  

    CMD="singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
            gatk --java-options \"-Xmx4g -Xms4g\" GenomicsDBImport \
                -V ./samples.list \
                --genomicsdb-workspace-path $OUT_DB \
                -L $INTERVAL \
                --reader-threads 4 \
                --batch-size 13"
      
    JOB_QSUB=$QSUB" -N $JID -o $JID.log -pe mpi 12"

    echo $CMD | $JOB_QSUB

done

#keep repeating until no failed jobs id'd







#genotype
for INTERVAL in $(cat ./contigs.list); do

    ############ genotype
    JID="geno_"$INTERVAL
    
    IN_DB=$(pwd)/db/$INTERVAL
    IN_REFERENCE=$HAE_GENOME
    OUT_VCF=$(pwd)/geno/$INTERVAL.vcf

    CMD="singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
            gatk GenotypeGVCFs \
                -R $IN_REFERENCE \
                -V gendb://$IN_DB \
                -new-qual \
                -O $OUT_VCF"
      
    JOB_QSUB=$QSUB" -N $JID -o logs/$JID.log -pe mpi 1"

    echo $CMD | $JOB_QSUB

    sleep 0.1
done


#check for failed jobs due to import issues
grep -l  "A USER ERROR has occurred: GenomicsDB workspace" logs/geno_*.log \
    | sed s'/logs\/geno_\(.*\).log/\1/gi' >failed_import.list


#continue to check for failed imports...then check for failed genotypes.
for INTERVAL in $(cat ./contigs.list); do

    if [ ! -f ./geno/$INTERVAL.vcf.idx ]; then
        echo $INTERVAL >>failed.list
    fi 

done

grep -L "Traversal complete" logs/geno_*.log | sed s'/logs\/geno_\(.*\).log/\1/gi' >>failed.list
grep -i -e kill -e die -e error logs/geno_*.log | sed s'/logs\/geno_\(.*\).log/\1/gi' >>failed.list
cut -f1 -d":" failed.list | sort -u >tmp
mv tmp failed.list

#so i re-ran several times and it looks like a
#re-run
for INTERVAL in $(cat ./failed.list); do

    ############ genotype
    JID="geno_"$INTERVAL
    
    IN_DB=$(pwd)/db/$INTERVAL
    IN_REFERENCE=$HAE_GENOME
    OUT_VCF=$(pwd)/geno/$INTERVAL.vcf

    if [ -f $OUT_VCF ]; then
        rm $OUT_VCF
    fi 

    if [ -f $OUT_VCF.idx ]; then
        rm $OUT_VCF.idx
    fi 

    if [ -f logs/$JID.log ]; then
        rm logs/$JID.log
    fi

    CMD="singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
            gatk GenotypeGVCFs \
                -R $IN_REFERENCE \
                -V gendb://$IN_DB \
                -new-qual \
                -O $OUT_VCF"
      
    JOB_QSUB=$QSUB" -N $JID -o logs/$JID.log -pe mpi 12"

    echo $CMD | $JOB_QSUB

done

#clean up import logs
for LOG in $(ls logs/import_*.log); do
    tar -rf gdbimport_logs.tar $LOG
    rm $LOG
done
gzip gdbimport_logs.tar

#clean up genotype logs
for LOG in $(ls logs/geno_*.log); do
    tar -rf genotype_logs.tar $LOG
    rm $LOG
done
gzip genotype_logs.tar







#merge all contig vcf files into a single cohort file
ls ./geno/*.vcf >vcf.list

#moving to titan
singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    split -n l/600 \
        --numeric-suffixes=1 \
        --additional-suffix .list \
        vcf.list \
        merge-1_


for CHUNK in $(seq -w 1 600); do
    CMD="singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
        gatk --java-options "-Xmx4g" MergeVcfs \
            --MAX_RECORDS_IN_RAM 500000 \
            -I merge-1_$CHUNK.list \
            -O merge-1_$CHUNK.vcf"

    echo $CMD | $QSUB -N merge-1_$CHUNK -o merge-1_$CHUNK.log -pe mpi 12
done



#final merge (on titan)
ssh titan
cd /master/nplatt/schisto_hybridization/results/call_whole_genome_snps
source activate snp_calling; source ../../scripts/set_env.sh

ls merge-1_*.vcf >vcf.list

singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    split -n l/10 \
        --numeric-suffixes=1 \
        --additional-suffix .list \
        vcf.list \
        merge-2_

for CHUNK in $(seq -w 1 10); do
    singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
        gatk --java-options "-Xmx8g" MergeVcfs \
            --MAX_RECORDS_IN_RAM 500000 \
            -I merge-2_$CHUNK.list \
            -O merge-2_$CHUNK.vcf \
            >merge-2_"$CHUNK".log 2>&1 &
done

#final merge
ls merge-2_*.vcf >vcf.list

singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    gatk --java-options "-Xmx50g" MergeVcfs \
        -I vcf.list \
        -O tmp.vcf

#sort
singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    gatk --java-options "-Xmx50g" SortVcf \
        -I tmp.vcf \
        -O sh_sb_wg_variants_raw.vcf


#clean dir
mv merge-*.log logs
rm merge-* tmp.vcf

################################################################################
#filter snps
#select snps and filter
singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    gatk SelectVariants \
        -V sh_sb_wg_variants_raw.vcf \
        -select-type SNP \
        -O wg_raw_snps.vcf \
        -R $HAE_GENOME \
        >select_snps.log 2>&1 &

singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    gatk SelectVariants \
        -V sh_sb_wg_variants_raw.vcf \
        -select-type INDEL \
        -O wg_raw_indels.vcf \
        -R $HAE_GENOME \
        >select_indels.log 2>&1 &

############
wait
############

singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    gatk VariantFiltration \
        -R $HAE_GENOME \
        -V wg_raw_snps.vcf \
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
        -O wg_soft-filtered_snps.vcf \
        >soft_filter_snps.log 2>&1 &


singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    gatk VariantFiltration \
        -R $HAE_GENOME \
        -V wg_raw_indels.vcf \
        --filter-name "QD_lt_2,indel" \
        --filter-expression "QD < 2.0" \
        --filter-name "FS_gt_200,indel" \
        --filter-expression "FS > 200.0" \
        --filter-name "ReadPosRankSum_lt_-20,indel" \
        --filter-expression "ReadPosRankSum < -20.0" \
        -O wg_soft-filtered_indels.vcf \
        >soft_filter_indels.log 2>&1 &

############
wait
############

#and merge them back together.
ls wg_soft-filtered_snps.vcf \
    wg_soft-filtered_indels.vcf \
    >vcf.list

wait 

singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    gatk MergeVcfs \
        -I vcf.list \
        -O wg_soft-filtered.vcf \
        -R $HAE_GENOME

#then hard filter
vcftools \
    --remove-filtered-all \
    --vcf wg_soft-filtered.vcf \
    --recode \
    --recode-INFO-all \
    --stdout \
    >shae-sbov_whole-genome_filtered_geno-qual.vcf
#After filtering, kept 13 out of 13 Individuals
#After filtering, kept 10,684,840 out of a possible 15,828,532 Sites

#index
${ENVIRONMENTS["SINGULARITY"]} \
     gatk IndexFeatureFile \
        -F shae-sbov_whole-genome_filtered_geno-qual.vcf




