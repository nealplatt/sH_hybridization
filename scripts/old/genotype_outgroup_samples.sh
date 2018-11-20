#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh

cd $GENO_ALL_DIR

mkdir $GENO_ALL_DIR/final_hc

########################
#VERY IMPORTANT TO NOTE.  
#ERR119622 and ERR119623 are long insert libraries to highly fragmented genomes
# as a result the paired reads do not map very well.  To accomdate this I am 
# re-running sections of the process BAM script to treat all reads as SE.

# <script here>

#CALL SNPS AT REF_PANEL LOCI IN EACH OF THE SAMPLES (SH AND OUTGROUP)
for SAMPLE in "${SAMPLES[@]}"; do

    SAMPLE=$(basename $BAM _processed.bam)

    echo"" >$SCRIPTS_DIR/$SAMPLE"_hc.sh"

    ############ HC
    JID=$SAMPLE"_hc"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD=""

    IN_BAM=$MAP_DIR/$SAMPLE"_processed.bam"
    OUT_VCF=$GENO_ALL_DIR/final_hc/$SAMPLE.vcf
    IN_REFERENCE=$HAE_GENOME
    IN_PANEL=$RESULTS_DIR/04-build_refpanel/sH_ref_panel_snps.list

    
    CMD="${ENVIRONMENTS[$ENV]} \
        gatk HaplotypeCaller \
            -I $BAM \
            -O $OUT_VCF \
            -R $IN_REFERENCE \
            -L $IN_PANEL \
            -ERC GVCF"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SCRIPTS_DIR/$SAMPLE"_hc.sh"
done

#WAIT UNTIL ALL JOBS ARE FINISHED.

#create list of vcf files from hc
ls $GENO_ALL_DIR/final_hc/*.vcf >$GENO_ALL_DIR/samples.list
cat $RESULTS_DIR/04-build_refpanel/sH_ref_panel_snps.list \
    | cut -f1 -d":" \
    | sort \
    | uniq \
    >$GENO_ALL_DIR/contigs.list

#make db dir
mkdir $GENO_ALL_DIR/final_db

#run GDIMPORT for each contig
for CONTIG in $(cat contigs.list); do

    echo"" >$SCRIPTS_DIR/$CONTIG"_import.sh"

    ############ GDBIMPORT
    JID=$CONTIG"_import"
    LOG=$LOGS_DIR/$JID".log"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD=""

    SAMPLES_LIST=$GENO_ALL_DIR/samples.list
    OUT_DB=$GENO_ALL_DIR/final_db/$CONTIG
    INTERVAL=$CONTIG
    
    CMD="${ENVIRONMENTS[$ENV]} \
            gatk --java-options \"-Xmx4g -Xms4g\" GenomicsDBImport \
                -V $SAMPLES_LIST \
                --genomicsdb-workspace-path $OUT_DB \
                -L $INTERVAL \
                --reader-threads $THREADS \
                --batch-size 24"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    #limit num running jobs to
    NUM_JOBS=$(qstat | grep nplatt | wc -l)

    while [ $NUM_JOBS -gt 330 ]; do
        echo "."
        sleep 1
        NUM_JOBS=$(qstat | grep nplatt | wc -l)
    done

    echo $CMD | $JOB_QSUB
    echo $CMD >>$SCRIPTS_DIR/$CONTIG"_import.sh"
done

#then genotype each sample
mkdir $GENO_ALL_DIR/final_genotype

for CONTIG in $(cat contigs.list); do

    ############ genotype
    JID=$CONTIG"_genotype"
    LOG=$LOGS_DIR/$JID".log"
    SCRIPT=$SCRIPTS_DIR/$JID".sh"
    THREADS=12
    ENV="SINGULARITY"     
    HOLD=""

    IN_DB=$GENO_ALL_DIR/final_db/$CONTIG
    IN_REFERENCE=$HAE_GENOME
    OUT_VCF=$GENO_ALL_DIR/final_genotype/$CONTIG.vcf

    CMD="${ENVIRONMENTS[$ENV]} \
            gatk GenotypeGVCFs \
                -R $IN_REFERENCE \
                -V gendb://$IN_DB \
                -new-qual \
                -O $OUT_VCF"
      
    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi $THREADS $HOLD"

    #limit num running jobs
    NUM_JOBS=$(qstat | grep nplatt | wc -l)

    while [ $NUM_JOBS -gt 400 ]; do
        echo "."
        sleep 1
        NUM_JOBS=$(qstat | grep nplatt | wc -l)
    done

    echo $CMD | $JOB_QSUB
    echo $CMD >$SCRIPT
done

#merge all contig vcf files into a single cohort file
ls $GENO_ALL_DIR/final_genotype/*.vcf >$GENO_ALL_DIR/vcf.list

${ENVIRONMENTS['SINGULARITY']} \
    split -n l/30 \
        --numeric-suffixes=1 \
        --additional-suffix .list \
        $GENO_ALL_DIR/vcf.list \
        $GENO_ALL_DIR/merge-01_


for CHUNK in $(seq -w 1 30); do
    CMD="${ENVIRONMENTS["SINGULARITY"]} \
        gatk --java-options "-Xmx4g" MergeVcfs \
            --MAX_RECORDS_IN_RAM 500000 \
            -I merge-01_$CHUNK.list \
            -O merge-01_$CHUNK.vcf"

    echo $CMD | $QSUB -N merge_$CHUNK -o merge_$CHUNK.log -pe mpi 12
done

ls $GENO_ALL_DIR/merge-01_*.vcf >$GENO_ALL_DIR/vcf.list
rm $GENO_ALL_DIR/merge-01_*.list

${ENVIRONMENTS["SINGULARITY"]} \
    gatk --java-options "-Xmx4g" MergeVcfs \
        --MAX_RECORDS_IN_RAM 500000 \
        -I $GENO_ALL_DIR/vcf.list \
        -O $GENO_ALL_DIR/tmp.vcf

${ENVIRONMENTS["SINGULARITY"]} \
    gatk --java-options "-Xmx8g" SortVcf \
        --MAX_RECORDS_IN_RAM 500000 \
        -I $GENO_ALL_DIR/tmp.vcf \
        -O $GENO_ALL_DIR/all_schisto_cohort_raw.vcf


#clean dir
rm $GENO_ALL_DIR/merge* $GENO_ALL_DIR/vcf.list $GENO_ALL_DIR/tmp.vcf*

#-------------------------------------------------------------------------------

mkdir final_filter
cd final_filter
#start filtering process on SNPS only
${ENVIRONMENTS["SINGULARITY"]} \
    gatk SelectVariants \
        -V $GENO_ALL_DIR/all_schisto_cohort_raw.vcf \
        -select-type SNP \
        -O all_schisto_cohort_raw_snps.vcf \
        -R $HAE_GENOME
#141703 of 151876 loci remaining

#calculate the percent missing per individual (we are going to save those with
# at least 66% of called snps
vcftools \
    --vcf all_schisto_cohort_raw_snps.vcf \
    --missing-indv \
    --stdout \
    >indiv_missing_table.tsv

cat indiv_missing_table.tsv \
    | awk '$5 >0.33 {{print $1}}' \
    | sed 1d \
    >remove_indiv.list
#have
#ERR103048 bovis
#ERR310937 curassoni
#ERR084970 haematobium
#ERR037800 haematobium
#SRR433865 haematobium
#ERR103051 bovis
#ERR119612 guineensis
#ERR119613 intercalatum
#ERR539855 mattheei
#ERR539857 mattheei
#ERR310940 margrebowiei

#lost
#ERR119622 bovis
#ERR119623 curassoni
#ERR539850 guineensis
#ERR539851 mattheei
#ERR539852 guineensis
#ERR539853 bovis
#ERR539854 intercalatum 
#ERR539856 intercalatum 
#Sh.TZ_PEM0075.1
#Sm.BR_0447.1
#Sm.BR_1278.1
#Sm.BR_2039.1


vcftools \
    --vcf all_schisto_cohort_raw_snps.vcf \
    --remove remove_indiv.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >all_schisto_cohort_raw_snps_66p.vcf
#After filtering, kept 105 out of 117 Individuals
#After filtering, kept 141703 out of a possible 141703 Sites

#soft filtered
${ENVIRONMENTS["SINGULARITY"]} \
    gatk VariantFiltration \
    -R $HAE_GENOME \
    -V all_schisto_cohort_raw_snps_66p.vcf \
    --filter-name "snp_QD_lt_5" \
    --filter-expression "QD < 5.0" \
    --filter-name "snp_FS_gt_55" \
    --filter-expression "FS > 55.0" \
    --filter-name "snp_MQ_lt_40" \
    --filter-expression "MQ < 40.0 || MQ = \"nan\"" \
    --filter-name "snp_MQRankSum_lt_-12.5" \
    --filter-expression "MQRankSum < -12.5" \
    --filter-name "snp_ReadPosRankSum_lt_-8" \
    --filter-expression "ReadPosRankSum < -8.0" \
    --filter-name "snp_SQR_gt_3" \
    --filter-expression "SOR > 3.0"  \
    -O all_schisto_cohort_soft.vcf

#hard filter 
vcftools \
    --remove-filtered-all \
    --min-alleles 2 \
    --max-alleles 4 \
    --maf 0.05 \
    --vcf all_schisto_cohort_soft.vcf \
    --recode \
    --recode-INFO-all \
    --stdout \
    >all_schisto_cohort_hard.vcf
#After filtering, kept 105 out of 105 Individuals
#After filtering, kept 129322 out of a possible 141703 Sites

#only want 2 missing since that means atleast 1 matthei or marg has call
vcftools \
    --max-missing-count 2 \
    --vcf all_schisto_cohort_hard.vcf \
    --missing-site \
    --indv ERR103048 \
    --indv ERR310937 \
    --indv ERR084970 \
    --indv ERR037800 \
    --indv SRR433865 \
    --indv ERR103051 \
    --indv ERR119612 \
    --indv ERR119613 \
    --indv ERR539855 \
    --indv ERR539857 \
    --indv ERR310940 \
    --stdout \
    | awk '{print $1"\t"$2}' \
    | sed 1d \
    >outgroup_keep_site.list

vcftools \
    --vcf all_schisto_cohort_hard.vcf \
    --positions outgroup_keep_site.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >all_schisto_cohort_9of11out.vcf
#After filtering, kept 105 out of 105 Individuals
#After filtering, kept 94559 out of a possible 129322 Sites

#now get snps that are single locus in mansoni genome.
#lift over to smansoni coordinates

#annotate each snp with a name - will be used to cull out multi pos later.
bcftools annotate \
    --set-id +'%CHROM\:%POS' \
    all_schisto_cohort_9of11out.vcf \
    >all_schisto_cohort_annotated.vcf

#convert to bed
vcf2bed \
    --do-not-sort \
    --max-mem=2G \
    <all_schisto_cohort_annotated.vcf \
    | awk '{print $1"\t"$2"\t"$3"\t"$4}' \
    >all_schisto_cohort_annotated.bed
    #95,165 remaining snps

#lift to sman coordinates
${ENVIRONMENTS["SINGULARITY"]} \
    /usr/software/progressiveCactus/submodules/hal/bin/halLiftover  \
        $RESULTS_DIR/06-wga/schMan7_vs_schMan1.hal \
        schHae_v1 \
        all_schisto_cohort_annotated.bed \
        schMan_v7 \
        all_schisto_cohort_sman.bed
        #109,507 remaining snps

#now lift the vcf to the sman coords (this will remove snp loci not aligned to sman)
python $WORK_DIR/scripts/lift_over_vcf.py \
    all_schisto_cohort_sman.bed \
    all_schisto_cohort_annotated.vcf \
    all_schisto_cohort_sman.vcf
    #81,394 remaining snps

#find loci that map to only a single region
awk '{print $4}' all_schisto_cohort_sman.bed \
    | sort \
    | uniq -c \
    | awk '$1==2 {print $2}' \
    >all_schisto_cohort_sman_multipos_remove.list
    #6,594 multiloci snps to remove


#now look for multiple snps aligned to a single locus
cut -f1,2,3 all_schisto_cohort_sman.bed \
    | sort \
    | uniq -c \
    | awk '$1==2 {print $2"\t"$3}' >all_schisto_cohort_sman_multipos_remove.pos
    #1,590 loci annotated wtih multiple snps

#and put them in a vcf
vcftools \
    --vcf all_schisto_cohort_sman.vcf \
    --exclude-positions all_schisto_cohort_sman_multipos_remove.pos \
    --exclude all_schisto_cohort_sman_multipos_remove.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >all_schisto_cohort_smansingle.vcf
    #After filtering, kept 105 out of 105 Individuals
    #After filtering, kept 74783 out of a possible 81394 Sites

#fix header and sort
grep -v "contig=<ID=" all_schisto_cohort_smansingle.vcf \
    >all_schisto_cohort_headerless.vcf

${ENVIRONMENTS["SINGULARITY"]} \
    gatk SelectVariants \
        -R $MAN_GENOME \
        -V all_schisto_cohort_headerless.vcf \
        -O all_schisto_cohort_headerfixed.vcf

${ENVIRONMENTS["SINGULARITY"]} \
    gatk SortVcf \
        -I all_schisto_cohort_headerfixed.vcf \
        -O ../all_schisto_cohort_filtered.vcf
        #74,783 in final (all chr) dataset
#-------------------------------------------------------------------------------

#get autosomal
vcftools \
    --vcf ../all_schisto_cohort_filtered.vcf \
    --chr SM_V7_1 \
    --chr SM_V7_2 \
    --chr SM_V7_3 \
    --chr SM_V7_4 \
    --chr SM_V7_5 \
    --chr SM_V7_6 \
    --chr SM_V7_7 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >autosomes_schisto_cohort_filtered.vcf
    #58,529 autosomal snps

#LD PRUNE
plink \
    --vcf autosomes_schisto_cohort_filtered.vcf \
    --allow-extra-chr \
    --indep-pairwise 25 5 0.2\
    --out autosomes_schisto_cohort_LD-25-5-2

vcftools \
    --vcf autosomes_schisto_cohort_filtered.vcf \
    --snps autosomes_schisto_cohort_LD-25-5-2.prune.in \
    --recode \
    --recode-INFO-all \
    --stdout \
    >../autosomes_schisto_cohort_filtered_LD.vcf
    #After filtering, kept 105 out of 105 Individuals
    #After filtering, kept 6571 out of a possible 58529 Sites

cd ..

grep "#" autosomes_schisto_cohort_filtered_LD.vcf  \
    | tail -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >autosomes_schisto_cohort_filtered_LD.samples

chmod a-w all_schisto_cohort_filtered.vcf
chmod a-w autosomes_schisto_cohort_filtered_LD.vcf

####################################################################################################

#FST between bovis and tanzania
cat $GENO_ALL_DIR/admixture_samples.list | grep "Sh.TZ" >$GENO_ALL_DIR/zanzibar.list                                       
echo "ERR103048" >$GENO_ALL_DIR/bovis.list
echo -e "ERR037800\nERR084970\nSRR433865" >$GENO_ALL_DIR/sra_haem.list


vcftools \
    --vcf all_schisto_cohort_raw_snps_75p_hard_annotated_single_smancoords_header_sorted_LD-25-5-2.vcf \
    --weir-fst-pop zanzibar.list \
    --weir-fst-pop bovis.list \
    --stdout \
    | awk '$3==1 {print $1"\t"$2}' \
        >bovis_zanzibar_fixed_snps.list                
#574 snps fst=1


vcftools \
    --vcf all_schisto_cohort_raw_snps_75p_hard_annotated_single_smancoords_header_sorted_LD-25-5-2.vcf \
    --weir-fst-pop sra_haem.list \
    --weir-fst-pop bovis.list \
    --stdout \
    | awk '$3==1 {print $1"\t"$2}' \
        >bovis_sra_haem_fixed_snps.list  
#2249 snps fst=1

#now get only those snps that are fixed between bovis/zanzibar
vcftools \
    --vcf all_schisto_cohort_raw_snps_75p_hard_annotated_single_smancoords_header_sorted_LD-25-5-2.vcf \
    --positions bovis_zanzibar_fixed_snps.list   \
    --stdout \
    --recode \
    >bovis_zanzibar_fixed_snps_all.vcf      

#samething with haem_sra
vcftools \
    --vcf all_schisto_cohort_raw_snps_75p_hard_annotated_single_smancoords_header_sorted_LD-25-5-2.vcf \
    --positions bovis_sra_haem_fixed_snps.list   \
    --stdout \
    --recode \
    >bovis_haem-sra_fixed_snps_all.vcf   

#convert this list to table for hindex/introgress





#mitochondria - pull from all filtered sites (since mito in LD)
vcftools \
    --vcf $GENO_ALL_DIR/all_schisto_smancoord_filtered.vcf \
    --chr SM_V7_MITO \
    --recode \
    --recode-INFO-all \
    --stdout \
    >$GENO_ALL_DIR/all_schisto_smancoord_filtered_mito.vcf
#After filtering, kept 94 out of 94 Individuals
#After filtering, kept 215 out of a possible 96071 Sites

 ../../scripts/vcf_to_diploid_fasta.py all_schisto_smancoord_filtered_mito.vcf all_schisto_smancoord_filtered_mito.fas
#still missing mito data for a sig number of NE samples.  



vcftools \
    --vcf $GENO_ALL_DIR/all_schisto_smancoord_filtered.vcf \
    --chr SM_V7_1 \
    --chr SM_V7_2 \
    --chr SM_V7_3 \
    --chr SM_V7_4 \
    --chr SM_V7_5 \
    --chr SM_V7_6 \
    --chr SM_V7_7


#**ERR037800 - haem - SAMEA980563 - unk
#**ERR084970 - haem - SAMEA1034756 - unk
#**ERR103048 - bovis - SAMEA1034751 - 	Iranga, Tanzania - ERS074976
#*ERR119622 - bovis - SAMEA1463529 - 	Iranga, Tanzania - ERS094852
#*ERR119623 - curassoni - SAMEA1034756 - Dakar, Senegal - unk
#ERR310937 - curassoni - ERS076740 - Dakar, Senegal
#**SRR433865 - haematobium - SAMN00794683 - Egypt

