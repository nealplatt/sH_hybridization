################################################################################
#  classify samples based on single gene hybrids
#
source activate hybrid_id

cd $RESULTS_DIR/07-single_gene_class

#get ITS and COI from each sample
    #where is it in mansoni
    #extract same region in samples
    #turn that region into a fasta
    
# compare to references
#-----------------------#
# REFERENCE SEQUENCES   #
#spec | accession | gene#
#-----------------------#
#bov	AY157212	COI
#bov	AJ519521	COI
#bov	FJ588862	ITS
#hae	GU257338	COI
#hae	GU257354	COI
#hae	JQ397385	COI
#hae	JQ397386	COI
#hae	JQ397401	ITS
#hae	JQ397406	ITS
#man	AY446106	COI
#man	JQ289619	COI
#man	JQ289753	ITS
#man	AF531314	ITS
#rod	AY446142	COI
#rod	AY446143	COI
#rod	AY446078	ITS
#rod	AF531312	ITS
#-----------------------

#1 blast against new sman assembly to ID coordinates (for extractinon)
#>JQ397386.1_schHae_cox
#AATGCCTGTATTGATAGGTGGATTTGGTAAATATTTTCTTCCGTTTTTTTTATATATAGATGATTTATTGTTACCTCGGTTGAATTCTTTTAGTTTATGATTAATGATTCCTTCATTTTTTTATATGGAGTTGAGTTTATATTATGGTTGTGGTGTAGGATGAACATTGTATCCTCCATTATCCATATCTGAGAATTCGGGTTTAGGTGTAGATTATTTAATGTTTTCTTTACATTTAGCGGGTGTATCTAGATTAGTTGGTTCCATTAATTTTATTTCTACGATTATTAGTCGTGTCAATTTTAAGACTTCTATAATAATATGATCATATTTGTTCACTTCTATCTTATTATTGCTTTCATTACCAGTTTTAGCAGCTGGTATTACTATGTTATTATTTGATCGTAAATTTGGTACTGCTTTTTTTGAGCCTATGGGTGGTGGTGATCCATTATTATTTCAGCACTTATTTTGATTTTTTGGTCATCCGGAGGTGTATGTTTTAATTTTACCTGGATTTGGAATAGTTAGTCATATATGTATGAGGATAAGTAATAATGATTCATCGTTTGGGTATTATGGATTAATTTGTGCTATGGCTTCGATAGTTTGTTTAGGAAGTGTAGTCTGAGCCCATCATATGTTTATGGTTGGTTTAGACTATTTGACTGCTATATTTTTTAGTTCAGTGACTATGATTATAGGGATTCCTACAGGTATAAAG-GTTTTTTCTTGATTATATATGCTTAAAAGCTGTGGGTCTCGTGTATGAGATCCTATAGTTTGATGATTGGTTGGTTTTATATTTTTATTTACGATAGGCGGTGTTACTGGTATAGCTTTATCAGCTTCTTCATTAGATATATTATTTCATGATACTTGATTTGTTGTTGCTCATTTTCATTATGTTCTTTCTTTAGGTTCTTATAGAAGTGTAGTAATAATGTTGTTATGGT
#>JQ397406.1_schHae_its
#TTATCATAACCCAAAATATATAATGATGCATGCACCTGGCTTCTTGCTGGACTGTATGTACCCTGGCTTGGTGGTTATTACCCTAGGCTTCAGTGGTTGATATTTTCTTGACCGGGGTACCTAGCCTGTCGTATGCCCTGATGGTGTTCTCGTAACTTTCGGGTTGCCTGATCTGCCAAGGGCGATGGGACAGTGCATGACGCTATTGTTGTGTGCTAGGTTCAAAGAGAATTGTATGCTATATGCATGCAAATCCGCCCCGTTATTGTTCCTATTTCAAACTTTTACACTGTTGAAGCGATCCGGTTTGGCTTGCCATTCACGGGTTTGCTGCCTGGCATGCACCTGGCTTCGTGCTGGACTGCATGTACGCTGGCTTAGCGGTAAATATCCTAGGCTGCAGCGTTAACCATTAGTTCTATGCATTTGGGAAACCAATGTATGGGATTATTGGCGTACAACTTTGAGCGGTGGATCACTCGGCTCGTGTGTCGATGAAGAGTGCAGCCAACTGTGTGAATTAATGTGAACTGCATACTGCTTTGAACATCGACATCTTGAACGCATATTGCGGCTACGGGATATCCTGTGGCCACGCCTGTCCGAGGGTCGGCTTTTCATCTATCACGGCGCACATTGAGTCGTGGATTGGGCGAGTGCCTGCCGGCGTGTATACCCGCATATCAACGCGGGTTGCTGGTCGAAGGCTCCGTCCTAATAGTCCGGCCACAGCCTAGTCCGGTCTAGATGACTTGATCGAGATGCTGCGGTGGGTTGTGCTCGAGTCGTGGCTTAATGACATTATACGCGCTCGGGAAGAATCGCACCTATCGTACGCTACGTTGGTCACTTGATCTTGTCTCTATGGTTCGGTCTACGGTTTGTACCGATGGTGTGTGTAATACGCACGAATTGTATAATTGACCC


echo ">JQ397386.1_schHae_cox" >schHae_ref_cox.fas
echo "AATGCCTGTATTGATAGGTGGATTTGGTAAATATTTTCTTCCGTTTTTTTTATATATAGATGATTTATTGTTACCTCGGTTGAATTCTTTTAGTTTATGATTAATGATTCCTTCATTTTTTTATATGGAGTTGAGTTTATATTATGGTTGTGGTGTAGGATGAACATTGTATCCTCCATTATCCATATCTGAGAATTCGGGTTTAGGTGTAGATTATTTAATGTTTTCTTTACATTTAGCGGGTGTATCTAGATTAGTTGGTTCCATTAATTTTATTTCTACGATTATTAGTCGTGTCAATTTTAAGACTTCTATAATAATATGATCATATTTGTTCACTTCTATCTTATTATTGCTTTCATTACCAGTTTTAGCAGCTGGTATTACTATGTTATTATTTGATCGTAAATTTGGTACTGCTTTTTTTGAGCCTATGGGTGGTGGTGATCCATTATTATTTCAGCACTTATTTTGATTTTTTGGTCATCCGGAGGTGTATGTTTTAATTTTACCTGGATTTGGAATAGTTAGTCATATATGTATGAGGATAAGTAATAATGATTCATCGTTTGGGTATTATGGATTAATTTGTGCTATGGCTTCGATAGTTTGTTTAGGAAGTGTAGTCTGAGCCCATCATATGTTTATGGTTGGTTTAGACTATTTGACTGCTATATTTTTTAGTTCAGTGACTATGATTATAGGGATTCCTACAGGTATAAAG-GTTTTTTCTTGATTATATATGCTTAAAAGCTGTGGGTCTCGTGTATGAGATCCTATAGTTTGATGATTGGTTGGTTTTATATTTTTATTTACGATAGGCGGTGTTACTGGTATAGCTTTATCAGCTTCTTCATTAGATATATTATTTCATGATACTTGATTTGTTGTTGCTCATTTTCATTATGTTCTTTCTTTAGGTTCTTATAGAAGTGTAGTAATAATGTTGTTATGGT" >schHae_ref_cox.fas    
echo ">JQ397406.1_schHae_its" >schHae_ref_its.fas
echo "TTATCATAACCCAAAATATATAATGATGCATGCACCTGGCTTCTTGCTGGACTGTATGTACCCTGGCTTGGTGGTTATTACCCTAGGCTTCAGTGGTTGATATTTTCTTGACCGGGGTACCTAGCCTGTCGTATGCCCTGATGGTGTTCTCGTAACTTTCGGGTTGCCTGATCTGCCAAGGGCGATGGGACAGTGCATGACGCTATTGTTGTGTGCTAGGTTCAAAGAGAATTGTATGCTATATGCATGCAAATCCGCCCCGTTATTGTTCCTATTTCAAACTTTTACACTGTTGAAGCGATCCGGTTTGGCTTGCCATTCACGGGTTTGCTGCCTGGCATGCACCTGGCTTCGTGCTGGACTGCATGTACGCTGGCTTAGCGGTAAATATCCTAGGCTGCAGCGTTAACCATTAGTTCTATGCATTTGGGAAACCAATGTATGGGATTATTGGCGTACAACTTTGAGCGGTGGATCACTCGGCTCGTGTGTCGATGAAGAGTGCAGCCAACTGTGTGAATTAATGTGAACTGCATACTGCTTTGAACATCGACATCTTGAACGCATATTGCGGCTACGGGATATCCTGTGGCCACGCCTGTCCGAGGGTCGGCTTTTCATCTATCACGGCGCACATTGAGTCGTGGATTGGGCGAGTGCCTGCCGGCGTGTATACCCGCATATCAACGCGGGTTGCTGGTCGAAGGCTCCGTCCTAATAGTCCGGCCACAGCCTAGTCCGGTCTAGATGACTTGATCGAGATGCTGCGGTGGGTTGTGCTCGAGTCGTGGCTTAATGACATTATACGCGCTCGGGAAGAATCGCACCTATCGTACGCTACGTTGGTCACTTGATCTTGTCTCTATGGTTCGGTCTACGGTTTGTACCGATGGTGTGTGTAATACGCACGAATTGTATAATTGACCC" >schHae_ref_its.fas

sed -i 's/-//g' schHae_ref_*.fas

makeblastdb -dbtype nucl -in $HAEMATOBIUM_GENOME

blastn \
    -subject $HAEMATOBIUM_GENOME \
    -query schHae_ref_cox.fas \
    -evalue 1e-10 \
    -out schHae_ref_cox_blastn.out \
    -outfmt 6 \
    -max_target_seqs 1 \
    -max_hsps 1

blastn \
    -subject $HAEMATOBIUM_GENOME \
    -query schHae_ref_its.fas \
    -evalue 1e-10 \
    -out schHae_ref_its_blastn.out \
    -outfmt 6 \
    -max_target_seqs 1 \
    -max_hsps 1


#get coordinates for sman
awk '{print $2":"$9"-"$10}' schHae_ref_cox_blastn.out #AMPZ01026399.1:2239-3194
awk '{print $2":"$9"-"$10}' schHae_ref_its_blastn.out #KL250502.1:40729-41391

#subsetting genome and reads to overlap with ITS and cox to speed up processeing
samtools faidx $HAEMATOBIUM_GENOME AMPZ01026399.1 KL250502.1 >schHae_ref_contigs.fa
${ENVIRONMENTS["SINGULARITY"]} gatk CreateSequenceDictionary -R schHae_ref_contigs.fa
samtools faidx schHae_ref_contigs.fa


        JID=$SAMPLE"_"$READ"_map"
        LOG=$LOGS_DIR/$JID".log"
        THREADS=12 
        ENV="CONDA_SNPS"
        HOLD="-hold_jid "$SAMPLE"_filter"
    
        SAI=$MAP_DIR/$SAMPLE"_"$READ".sai"
        FQ=$FILTER_DIR/$SAMPLE"_filtered_paired_"$READ".fastq.gz"
        GENOME=${GENOMES[$SAMPLE]}

        JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

        CMD="${ENVIRONMENTS[$ENV]} bwa aln -t 12 -f $SAI $GENOME $FQ"
            
        echo $CMD | $JOB_QSUB
        echo $CMD >>$SAMPLE"_process.sh"


#call snps on cox and its
mkdir cox-its_vcfs
for BAM in $(ls $RESULTS_DIR/01-processed_reads/sh_mapped_reads/*processed.bam); do

    SAMPLE=$(basename $BAM _processed.bam)

    JID=$SAMPLE"_cox-its_hc"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12

    #only call snps on target regions
    CMD="${ENVIRONMENTS["SINGULARITY"]} \
        gatk HaplotypeCaller \
            -I $BAM \
            -O cox-its_vcfs/$SAMPLE"_cox-its_bp.vcf" \
            -R $HAEMATOBIUM_GENOME  \
            -L AMPZ01026399.1:2239-3194 \
            -L KL250502.1:40729-41391 \
            -ERC BP_RESOLUTION"

    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"
            
    echo $CMD | $JOB_QSUB
done

#in the meantime call snps on the whole mitochondrial genome
mkdir whole_mito_vcfs
for BAM in $(ls $RESULTS_DIR/01-processed_reads/sh_mapped_reads/*processed.bam); do

    SAMPLE=$(basename $BAM _processed.bam)

    JID=$SAMPLE"_mito_hc"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12

    #only call snps on target regions
    CMD="${ENVIRONMENTS["SINGULARITY"]} \
        gatk HaplotypeCaller \
            -I $BAM \
            -O whole_mito_vcfs/$SAMPLE"_mito.g.vcf" \
            -R $HAEMATOBIUM_GENOME  \
            -L AMPZ01026399.1 \
            -ERC BP_RESOLUTION"

    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS"
            
    echo $CMD | $JOB_QSUB
done

ls cox-its_vcfs/*cox-its_bp.vcf >vcf.list
ls whole_mito_vcfs/*.vcf >mito_vcf.list


mkdir cox-its_hcdb

#run gdbimport (per contig)
${ENVIRONMENTS["SINGULARITY"]} \
    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
        -V mito_vcf.list \
        --genomicsdb-workspace-path mito_hcdb \
        -L AMPZ01026399.1 \
        --reader-threads 12 \
        --batch-size 12

${ENVIRONMENTS["SINGULARITY"]} \
    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
        -V vcf.list \
        --genomicsdb-workspace-path cox-its_hcdb/cox \
        -L AMPZ01026399.1 \
        --reader-threads 12 \
        --batch-size 12

${ENVIRONMENTS["SINGULARITY"]} \
    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
        -V vcf.list \
        --genomicsdb-workspace-path cox-its_hcdb/its \
        -L KL250502.1:40729-41391 \
        --reader-threads 12 \
        --batch-size 12


################################################################################
${ENVIRONMENTS["SINGULARITY"]} \
    gatk GenotypeGVCFs \
        -R $HAEMATOBIUM_GENOME  \
        -V gendb://cox-its_hcdb/its \
        -new-qual \
        -O its.bp.vcf

${ENVIRONMENTS["SINGULARITY"]} \
    gatk GenotypeGVCFs \
        -R $HAEMATOBIUM_GENOME  \
        -V gendb://cox-its_hcdb/cox \
        -new-qual \
        -O cox.bp.vcf

${ENVIRONMENTS["SINGULARITY"]} \
    gatk GenotypeGVCFs \
        -R $HAEMATOBIUM_GENOME  \
        -V gendb://mito_hcdb \
        -new-qual \
        -O mito.bp.vcf



vcftools \
    --vcf cox.bp.vcf \
    --recode \
    --remove-filtered-all \
    --remove-indels \
    --min-alleles 2 \
    --max-alleles 2 \
    --max-missing 0.51 \
    --maf 0.05 \
    --stdout >cox_filtered.vcf
#35 snps


vcftools \
    --vcf mito.bp.vcf \
    --recode \
    --remove-filtered-all \
    --remove-indels \
    --min-alleles 2 \
    --max-alleles 2 \
    --max-missing 0.51 \
    --maf 0.05 \
    --stdout >mito_filtered.vcf

vcftools \
    --vcf its.bp.vcf \
    --recode \
    --remove-filtered-all \
    --remove-indels \
    --min-alleles 2 \
    --max-alleles 2 \
    --max-missing 0.51 \
    --maf 0.05 \
    --stdout >its_filtered.vcf

python ../04-popstats/vcf_to_diploid_fasta.py mito_filtered.vcf mito_filtered.fas
python ../04-popstats/vcf_to_diploid_fasta.py cox_filtered.vcf cox_filtered.fas
python ../04-popstats/vcf_to_diploid_fasta.py its_filtered.vcf its_filtered.fas


vcftools \
    --vcf its.g.vcf \
    --chr KL250502.1 \
    --from-bp 40712 \
    --to-bp 41391 \
    --remove-filtered-all \
    --min-alleles 2 \
    --max-alleles 2 \
    --max-missing 0.51 \
    --maf 0.05 \
    --recode \
    --stdout >its_target.vcf
#14 snps


vcftools \
    --vcf its_filtered.vcf \
    --missing-indv \
    --stdout >its.missing_indiv




        #remove individuals missing more than 25% of the data    
        cat {output.SEX_MISSING_TABLE} | awk '$5 >0.25 {{print $1}}' >{output.SEX_REMOVE_LIST}

        vcftools \
            --vcf {output.SEX_FILT_1} \
            --remove {output.SEX_REMOVE_LIST} \
            --recode \
            --recode-INFO-all \
            --stdout >{output.SEX_FILT_2}





vcftools \
    --bed sman_its_gene_coords.bed \
    --vcf ../04-popstats/sH_filtered_auto-variants_sManLift.vcf \
    --recode \
    --stdout >its.vcf

vcftools \
    --bed sman_cox_gene_coords.bed \
    --vcf ../04-popstats/sH_filtered_mito-variants_sManLift.vcf \
    --recode \
    --stdout >cox.vcf


#convert vcf to fasta

#examine against reference panel
