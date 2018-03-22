source /master/nplatt/sH_hybridization/scripts/set-env.sh

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
$SINGULARITY bedtools complement -i all_mito_reads_merged.bed -g schHae_v1_mt.fa.fai >mito_to_exclude_no_baits.bed

#have to physically make a genome file 
#AMPZ01026399.1  17526
$SINGULARITY bedtools complement -i all_mito_reads_merged.bed -g schHae_v1_mt.fa.fai >mito_to_exclude_no_baits.bed

#now mask file and create necessary indecies
$SINGULARITY bedtools maskfasta -fi schHae_v1_mt.fa -fo schHae_v1_mt.masked.fa -bed mito_to_exclude_no_baits.bed
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


    java -jar ~/bin/gatk-3.8.0.jar \
        -T FastaAlternateReferenceMaker \
        -R schHae_v1_mt.masked.fa  \
        -o $INDIVIDUAL.mito.fas \
        -V mito_variants.vcf

    #change the sequence header name
    sed -i 's/>.*/>'$INDIVIDUAL'##AMPZ01026399.1/' $INDIVIDUAL.mito.fas

done

#create a fasta file with all mito sequences
cat *mito.fas >sh.mito.fas


#get samples form ncbi as outgroups

#bovis genome
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/HGI/Schistosoma_bovis_v1.0.fa.gz


#align

#exclude regions not in vcf
rm all_mito_reads.bed
for BAM in $(ls ../map_reads/*processed.bam); do
    $SINGULARITY bedtools bamtobed -i $BAM | grep AMPZ01026399.1 >>all_mito_reads.bed
done



 


