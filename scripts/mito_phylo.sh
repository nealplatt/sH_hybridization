source /master/nplatt/sH_hybridization/scripts/set-env.sh

mkdir mito_phylo

cd mito_phylo

singularity exec ../../snpCalling_v0.0.7.img \
    gatk SelectVariants \
        -R $REFERENCE \
        -V ../filter_cohort_vcf/sHaem_filtered.vcf  \
        -L AMPZ01026399.1 \
        -O mito_variants.vcf
#only 2 -- probably not worth doing unless re-filtering


for INDIVIDUAL in $(cat ../sample.list); do
    $SINGULARITY gatk SelectVariants \
        -R $REFERENCE \
        -V ../filter_cohort_vcf/sHaem_filtered.vcf \
        -sn $INDIVIDUAL \
        -O $INDIVIDUAL"_filtered.vcf"
done

#FILTER AND MERGE VARIANTS FROM MITOCHONDRIA
for SAMPLE in $(cat $SAMPLE_LIST); do
    # SELECT SNPS --------------------------------------------------------------
    JOB_NAME="mito_"$SAMPLE"_select_snps"
    THREADS=1
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND=""
    SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    IN_VCF="../haplotype_caller/"$SAMPLE"_final.g.vcf"
    OUT_VCF=$SAMPLE"_mito_raw-snps.g.vcf"
    
    CMD="$SINGULARITY gatk SelectVariants \
        -V $IN_VCF \
        -L AMPZ01026399.1 \
        -O $OUT_VCF \
        -R $REFERENCE"

    #DELETE $LOG $SCRIPT
    #SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"

    # FILTER SNPS --------------------------------------------------------------
    JOB_NAME="mito_"$SAMPLE"_filter_snps"
    THREADS=1
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid mito_"$SAMPLE"_select_snps"
    SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    IN_VCF=$OUT_VCF
    OUT_VCF=$SAMPLE"_mito_filtered-snps.g.vcf"

    CMD="$SINGULARITY gatk VariantFiltration \
        -R $REFERENCE \
        -V $IN_VCF \
        --filter-name "MIN_DP" \
        --filter-expression "DP > 100" \
        --filter-expression "QD < 2.0" \
        --filter-expression "MQ < 40.0" \
        --filter-expression "MQRankSum < -12.5" \
        --filter-expression "ReadPosRankSum" \
        --filter-expression "DP > 100" \
        --filter-expression "DP > 100" \


--filter-expression "'"DP > 100 && QD < 2.0 && FS > 60.0 && MQ < 40.0 && MQRankSum < -12.5 && ReadPosRankSum < -8.0"'" \
        --filter-name "'"recommended_snp_filter"'" \
        -O $OUT_VCF"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"

    # SELECT INDELS ------------------------------------------------------------
    JOB_NAME="mito_"$SAMPLE"_select_indels"
    THREADS=1
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND=""
    SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    IN_VCF="../haplotype_caller/"$SAMPLE"_final.g.vcf"
    OUT_VCF=$SAMPLE"_mito_raw-indels.g.vcf"
    
    CMD="$SINGULARITY gatk SelectVariants \
        -V $IN_VCF \
        -select-type INDEL \
        -L AMPZ01026399.1 \
        -O $OUT_VCF \
        -R $REFERENCE"

    #DELETE $LOG $SCRIPT
    #SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"


    # FILTER INDELS ------------------------------------------------------------
    JOB_NAME="mito_"$SAMPLE"_filter_indels"
    THREADS=1
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid mito_"$SAMPLE"_select_indels"
    SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    IN_VCF=$OUT_VCF
    OUT_VCF=$SAMPLE"_mito_filtered-indels.g.vcf"

    CMD="$SINGULARITY gatk VariantFiltration \
        -R $REFERENCE \
        -V $IN_VCF \
        --filter-expression "'"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"'" \
        --filter-name "'"recommended_indel_filter"'" \
        -O $OUT_VCF"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"


    # MERGE VARIANTS------------------------------------------------------------
    JOB_NAME="snp.merge_filtered_variants"
    THREADS=1
    LOG="$LOGS_DIR/$JOB_NAME.log" 
    DEPEND="-hold_jid mito_"$SAMPLE"_select_indels,mito_"$SAMPLE"_select_snps"
    SCRIPT="$SUB_SCRIPTS_DIR/$JOB_NAME.sh"

    JOB_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o $LOG $DEPEND"

    echo -e $SAMPLE"_mito_filtered-indels.g.vcf\n"$SAMPLE"_mito_filtered-indels.g.vcf" >$SAMPLE.list

    IN_LIST=$SAMPLE.list
    OUT_VCF=$SAMPLE"_mito_variants.vcf"

    CMD="$SINGULARITY gatk MergeVcfs -I $IN_LIST -O $OUT_VCF -R $REFERENCE"

    DELETE $LOG $SCRIPT
    SUBMIT "$CMD" "$SCRIPT" "$JOB_QSUB"
done





IN_RAW_VCF=../base_recalibration/r3_vcfs/cohort_raw_r3.vcf

SNP_VCF=cohort_raw_SNPs_r3.vcf
INDEL_VCF=cohort_raw_INDELs_r3.vcf
REFERENCE=/master/nplatt/sH_hybridization/data/genome/schHae_v1.fa

FILTERED_INDELS_VCF=cohort_filtered_INDELs_r3.vcf
FILTERED_SNPS_VCF=cohort_filtered_SNPs_r3.vcf

MERGED_VARIANTS_VCF=cohort_filtered_variants_r3.vcf


#get mito seq
samtools faidx ../../data/genome/schHae_v1.fa AMPZ01026399.1 >schHae_v1_mt.fa 
$SINGULARITY gatk CreateSequenceDictionary --REFERENCE schHae_v1_mt.fa
samtools faidx schHae_v1_mt.fa 

#find a way to exclude sites that weren't in the vcf file
$SINGULARITY bedtools sort -i all_mito_reads.bed >tmp
$SINGULARITY bedtools merge -i tmp >all_mito_reads_merged.bed 
#have to physically make a genome file 
#AMPZ01026399.1  17526
$SINGULARITY bedtools complement -i all_mito_reads_merged.bed -g schHae_v1_mt.fa.fai >exclude.bed

#mask
$SINGULARITY bedtools maskfasta -fi schHae_v1_mt.fa -fo schHae_v1_mt.masked.fa -bed exclude.bed

#now mask file and create necessary indecies
$SINGULARITY gatk CreateSequenceDictionary --REFERENCE schHae_v1_mt.masked.fa
samtools faidx schHae_v1_mt.masked.fa 

#get mitochondrial snps
$SINGULARITY gatk SelectVariants \
     -R $REFERENCE \
     -V $MERGED_VARIANTS_VCF  \
     -L AMPZ01026399.1 \
     -O mito_variants.vcf


for INDIVIDUAL in $(cat ../sample.list); do
    $SINGULARITY gatk SelectVariants \
        -R $REFERENCE \
        -V mito_variants.vcf \
        -sn $INDIVIDUAL \
        -O $INDIVIDUAL.mito.vcf
done
    #python vcf2fasta.py -v $INDIVIDUAL.mito.vcf -o $INDIVIDUAL.mito.fas -c AMPZ01026399.1

    java -jar ~/bin/gatk-3.8.0.jar \
       -T FastaAlternateReferenceMaker \
       -R schHae_v1_mt.masked.fa  \
       -o $INDIVIDUAL.mito.fas \
       -V mito_variants.vcf

    #change the sequence header name
    #sed -i 's/>.*/>'$INDIVIDUAL'##AMPZ01026399.1/' $INDIVIDUAL.mito.fas

done

#create a fasta file with all mito sequences
cat *mito.fas >sh.mito.fas


#get samples form ncbi as outgroups
#bovis
#mansoni
#curasoni
#japonicum
#haem

#bovis genome
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/HGI/Schistosoma_bovis_v1.0.fa.gz


#align

#exclude regions not in vcf
rm all_mito_reads.bed
for BAM in $(ls ../map_reads/*processed.bam); do
    $SINGULARITY bedtools bamtobed -i $BAM | grep AMPZ01026399.1 >>all_mito_reads.bed
done



 


