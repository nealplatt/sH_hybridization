#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# get_schHaem_invadolysin_seqs.sh - build a fasta file with the phased
#       Smp127030 (invadolysin) sequences.

# Uses a conda and singularity to manage the enviroment; needs relativley high
#       mem compute (used Titan - 125Gb)

#set up environment
source activate snp_calling
source /master/nplatt/schisto_hybridization/scripts/set_env.sh

cd $RESULTS_DIR

mkdir get_Smp127030_sequences

cd get_Smp127030_sequences

grep cds ../../data/Sm_v7.0.bed | grep Smp_127030 >Smp_127030_cds.bed


#lift the coordinates to haematobium
/opt/projects/singularity-2.4.2/bin/singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    /usr/software/progressiveCactus/submodules/hal/bin/halLiftover  \
        $RESULTS_DIR/wga/schMan7_vs_schMan1.hal \
        schMan_v7 \
        Smp_127030_cds.bed \
        schHae_v1 \
        Smp_127030_cds_schHae.bed
 
awk '{SUM+=$3-$2} END {print SUM/3}' Smp_127030_cds_schHae.bed
#check to see that things lifted over correctly

#get a fasta and mask for uprobed regions
samtools faidx ../../data/genome/schHae_v1.fa KL250964.1 >KL250964.1.fasta


bedtools intersect \
    -v \
    -a Smp_127030_cds_schHae.bed \
    -b ../../data/schHae_v1_probes.bed  \
    >Smp_127030_UNprobed-cds_schHae.bed

bedtools maskfasta \
    -fi KL250964.1.fasta \
    -bed Smp_127030_UNprobed-cds_schHae.bed \
    -fo KL250964.1_masked.fasta

#index sequences
samtools faidx KL250964.1_masked.fasta

/opt/projects/singularity-2.4.2/bin/singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    gatk CreateSequenceDictionary \
    -R KL250964.1_masked.fasta


#then fasta from reference of invadolysin snps
#...need to convert vcf from sman to shaem coords

#get a small region of the vcf file (wiht sman coordinates)
vcftools \
    --vcf ../beagle/auto_beagle_maf05.vcf \
    --chr SM_V7_4 \
    --from-bp 19000000 \
    --to-bp 21000000 \
    --recode \
    --stdout \
    >auto_beagle_maf05_chr4:19M-21M.vcf 

#now convert these coordinates to schHaem
python ../../scripts/schMan_to_schHae_vcf_coords.py \
    auto_beagle_maf05_chr4:19M-21M.vcf  \
    tmp.vcf 

#make the vcf compatible by fixing the header
grep -v "contig=<ID=" tmp.vcf \
    >headerless.vcf

#add contigs for sman to header
/opt/projects/singularity-2.4.2/bin/singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    gatk SelectVariants \
        -R $HAE_GENOME \
        -V headerless.vcf \
        -O header.vcf

#may need to be run on high mem nodes
/opt/projects/singularity-2.4.2/bin/singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    gatk SortVcf \
        -R $HAE_GENOME \
        -I header.vcf \
        -O Smp_127030_cds_schHae_coords.vcf

#clean up
rm header*
rm tmp.vcf
rm out.log
rm auto_beagle_maf05_chr4:19M-21M.vcf


#now I need phased haplotypes for each sample...to do this need to run beagle
CHR=SM_V7_4

#extract autosome specific vcf for the samples of interest
vcftools \
    --vcf ../build_snp_panel/auto.vcf \
    --chr $CHR \
    --recode \
    --stdout >$CHR".vcf"

sed -i 's/,assembly=schMan_v7.fa//gi' $CHR.vcf

grep -v "#" ../build_snp_panel/auto.vcf  \
    | grep -i sm_v7_4 \
        | awk '{printf "%s\t%s\t%.6f\t%s\n", $1, $3, $2/287000, $2}' \
            >$CHR".map"

beagle \
    gt=$CHR".vcf" \
    out=$CHR"_beagle_all_samples" \
    map=$CHR".map" \
    nthreads=4 \
    window=300 \
    overlap=30 \
    niterations=250


#now get the phased invadolysin haplotypes
vcftools \
    --vcf chr4.vcf \
    --bed Smp_127030_cds.bed \
    --recode \
    --recode-INFO-all \
    --stdout \
    >Smp_127030_phased.vcf


#now create a fake file that generates on "haplotype" per vcf
python ../../scripts/diploid_to_haploid_vcf.py \
    Smp_127030_phased.vcf \
    Smp_127030_phased_hap_A.vcf \
    Smp_127030_phased_hap_B.vcf

#convert each haplotype to a schMan invadolysin vcf file
for HAPLOTYPE in A B; do
    
    #in out files
    IN_VCF="Smp_127030_phased_hap_"$HAPLOTYPE".vcf"
    OUT_VCF="Smp_127030_hap_"$HAPLOTYPE"_schHae_coords.vcf"

    #lift coordinates
    python ../../scripts/schMan_to_schHae_vcf_coords.py $IN_VCF tmp.vcf 

    #strip header
    grep -v "contig=<ID=" tmp.vcf >headerless.vcf

    #add contigs for schHae to header
    ${ENVIRONMENTS["TITAN SINGULARITY"]} \
        gatk SelectVariants \
            -R $HAE_GENOME \
            -V headerless.vcf \
            -O header.vcf

    #sort to generate final vcf file
    ${ENVIRONMENTS["TITAN SINGULARITY"]} \
        gatk SortVcf \
            -R $HAE_GENOME \
            -I header.vcf \
            -O $OUT_VCF

    #clean up
    rm tmp.vcf 
    rm headerless.vcf 
    rm header.vcf*
    rm out.log

done


#use haplotypededed vcf file to make a fasta file
for SAMPLE in $(cat ../samples.list); do
    for HAPLOTYPE in A B; do

        #in out files
        IN_VCF="Smp_127030_hap_"$HAPLOTYPE"_schHae_coords.vcf"
        
        #get the vcf
        vcftools \
            --vcf $IN_VCF \
            --indv $SAMPLE \
            --recode \
            --recode-INFO-all \
            --stdout \
            >sample.vcf

        #get the fasta
        gatk \
           -T FastaAlternateReferenceMaker \
           -R KL250964.1_masked.fasta \
           -IUPAC $SAMPLE\
           -o sample.fasta \
           -V sample.vcf

        #clean up fasta for extraction
        sed -i "1c>KL250964.1" sample.fasta
        samtools faidx sample.fasta

        #extract each cds and them combine into one
        bedtools getfasta \
            -bed Smp_127030_cds_schHae.bed \
            -fi sample.fasta \
            -fo exons.fasta

        #combine into a single sequence
         grep -v "^>" exons.fasta \
            | awk 'BEGIN { ORS=""; print ">Sequence_name\n" } { print }' > cds.fasta

        #change the header
        sed -i "1c>$HAPLOTYPE.$SAMPLE" cds.fasta
        
        #add it to combined fasta file
        cat cds.fasta >>Smp_127030_haplotypes_cds.fasta
        echo >>Smp_127030_haplotypes_cds.fasta

        #clean up
        rm sample.fasta
        rm sample.vcf*
        rm exons.fasta*
        rm cds.fasta
        rm out.log
    done
done






















#for each sample make a fasta
for SAMPLE in $(cat ../samples.list); do
    
    #get the vcf
    vcftools \
        --vcf Smp_127030_cds_schHae_coords.vcf \
        --indv $SAMPLE \
        --chr KL250964.1 \
        --recode \
        --recode-INFO-all \
        --stdout \
        >$SAMPLE"_Smp_127030_cds.vcf"

    #get the fasta
    java -jar ../../scripts/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar  \
       -T FastaAlternateReferenceMaker \
       -R KL250964.1_masked.fasta \
       -IUPAC $SAMPLE\
       -o $SAMPLE"_Smp_127030_cds.fasta" \
       -V $SAMPLE"_Smp_127030_cds.vcf"

    #clean up fasta for extraction
    sed -i "1c>KL250964.1" $SAMPLE"_Smp_127030_cds.fasta"
    samtools faidx $SAMPLE"_Smp_127030_cds.fasta"

    #extract each cds and them combine into one
    bedtools getfasta \
        -bed Smp_127030_cds_schHae.bed \
        -fi $SAMPLE"_Smp_127030_cds.fasta" \
        -fo $SAMPLE"_Smp_127030_cds_indiv.fasta"

    #combine into a single sequence
     grep -v "^>" $SAMPLE"_Smp_127030_cds_indiv.fasta" \
        | awk 'BEGIN { ORS=""; print ">Sequence_name\n" } { print }' > $SAMPLE"_Smp_127030_combined_cds.fasta"

    #change the header
    sed -i "1c>$SAMPLE" $SAMPLE"_Smp_127030_combined_cds.fasta"
    
    #add it to combined fasta file
    cat $SAMPLE"_Smp_127030_combined_cds.fasta" >>Smp_127030_cds.fasta
    echo >>Smp_127030_cds.fasta

    #clean up
    rm $SAMPLE"_Smp_127030_cds.vcf"
    rm $SAMPLE"_Smp_127030_cds.vcf.idx"
    rm $SAMPLE"_Smp_127030_cds.fasta"
    rm $SAMPLE"_Smp_127030_cds.fasta.fai"
    rm $SAMPLE"_Smp_127030_cds_indiv.fasta"
    rm $SAMPLE"_Smp_127030_combined_cds.fasta"
    rm out.log

done





