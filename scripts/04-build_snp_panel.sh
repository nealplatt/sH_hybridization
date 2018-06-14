mkdir $RESULTS_DIR/build_snp_panel
cd $RESULTS_DIR/build_snp_panel

#ERR119622 bovis
#ERR103048 bovis
#ERR310937 curassoni
#ERR119623 curassoni
#ERR084970 haematobium
#ERR037800 haematobium
#SRR433865 haematobium
#ERR103051 bovis
#ERR119612 guineensis
#ERR119613 intercalatum
#ERR310940 margrebowiei
#ERR539850 guineensis
#ERR539851 mattheei
#ERR539852 guineensis
#ERR539853 bovis
#ERR539854 intercalatum 
#ERR539855 mattheei
#ERR539856 intercalatum 
#ERR539857 mattheei

ln -s $SNP_DIR/cohort_raw_bqsr-1.vcf
ln -s $SNP_DIR/cohort_raw_bqsr-1.vcf.idx

${ENVIRONMENTS["SINGULARITY"]} \
    gatk SelectVariants \
        -V cohort_raw_bqsr-1.vcf \
        -select-type SNP \
        -O cohort_snps.vcf \
        -R $HAE_GENOME
    #** 611,210 SNPs **#

#annotate snps names
bcftools annotate \
    --set-id +'%CHROM\:%POS' \
    cohort_snps.vcf \
    >cohort_snps_annotated.vcf

#filter loci missing genotype calls in 33% of samples
vcftools \
    --max-missing 0.33 \
    --vcf cohort_snps_annotated.vcf \
    --missing-site \
     --stdout \
    | awk '{print $1"\t"$2}' \
    | sed 1d \
    >keep_site.list

vcftools \
    --vcf cohort_snps_annotated.vcf \
    --positions keep_site.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_67perSite.vcf
#610,745 SNPs
#-------------------------

#filter individuals not genotyped at 20% of sites.
#   80% genotyping rate
vcftools \
    --vcf cohort_snps_67perSite.vcf \
    --missing-indv \
    --stdout \
    >indiv_missing_table.tsv

cat indiv_missing_table.tsv \
    | awk '$5 >0.20 {{print $1}}' \
    | sed 1d \
    >remove_indiv.list

vcftools \
    --vcf cohort_snps_67perSite.vcf \
    --remove remove_indiv.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_67perIndiv.vcf
#610,745 SNPs
#X Individuals (from X)
#-------------------------


#find only the snps where one margeboweri/mattehi is genotyped
vcftools \
    --max-missing-count 2 \
    --vcf cohort_snps_67perIndiv.vcf \
    --missing-site \
    --indv ERR310940 \
    --indv ERR539855 \
    --indv ERR539857 \
    --stdout \
    | awk '{print $1"\t"$2}' \
    | sed 1d \
    >outgroup_keep_site.list

vcftools \
    --vcf cohort_snps_67perIndiv.vcf \
    --positions outgroup_keep_site.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_oneOutGeno.vcf
    #After filtering, kept 534,853 out of a possible 610,745 Sites
#-------------------------

#only keep bi-allelic sites
bcftools view \
    -m2 \
    -M2 \
    -v snps \
    cohort_snps_oneOutGeno.vcf \
    >cohort_snps_oneOutGeno_biallelic.vcf
    #After filtering, kept 504,072 sites

#now get snps that are single locus in mansoni genome.
#lift over to smansoni coordinates

#convert to bed
vcf2bed \
    --do-not-sort \
    --max-mem=2G \
    <cohort_snps_oneOutGeno_biallelic.vcf \
    | awk '{print $1"\t"$2"\t"$3"\t"$4}' \
    >cohort_snps_oneOutGeno_biallelic.bed
    #504,072 remaining snps

#lift to sman coordinates
${ENVIRONMENTS["SINGULARITY"]} \
    /usr/software/progressiveCactus/submodules/hal/bin/halLiftover  \
        $RESULTS_DIR/wga/schMan7_vs_schMan1.hal \
        schHae_v1 \
        cohort_snps_oneOutGeno_biallelic.bed \
        schMan_v7 \
        cohort_snps_schMan.bed
        #579,260 remaining snps (increase is due to multilocus snps in mansoni genome)

#now lift the vcf to the sman coords (this will remove snp loci not aligned to sman)
python $WORK_DIR/scripts/lift_over_vcf.py \
    cohort_snps_schMan.bed \
    cohort_snps_oneOutGeno.vcf \
    cohort_snps_schMan.vcf
    #487,655 remaining snps

#find loci that map to only a single region
awk '{print $4}' cohort_snps_schMan.bed \
    | sort \
    | uniq -c \
    | awk '$1>1 {print $2}' \
    >cohort_snps_schMan_multipos_remove.list
    #45,392 multiloci snps to remove


#now look for multiple snps aligned to a single locus
cut -f1,2,3 cohort_snps_schMan.bed \
    | sort \
    | uniq -c \
    | awk '$1>1 {print $2"\t"$3}' >cohort_snps_schMan_multipos_remove.pos
    #4,287 loci annotated wtih multiple snps

#and put them in a vcf
vcftools \
    --vcf cohort_snps_schMan.vcf \
    --exclude-positions cohort_snps_schMan_multipos_remove.pos \
    --exclude cohort_snps_schMan_multipos_remove.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_schMan_singlePos.vcf
    #After filtering, kept 105 out of 105 Individuals
    #After filtering, kept 442,247 sites


#fix header and sort
grep -v "contig=<ID=" cohort_snps_schMan_singlePos.vcf \
    >cohort_snps_schMan_headerless.vcf

${ENVIRONMENTS["SINGULARITY"]} \
    gatk SelectVariants \
        -R $MAN_GENOME \
        -V cohort_snps_schMan_headerless.vcf \
        -O cohort_snps_schMan_header.vcf

#may need to be run on titan
${ENVIRONMENTS["SINGULARITY"]} \
    gatk SortVcf \
        -R $MAN_GENOME \
        -I cohort_snps_schMan_header.vcf \
        -O cohort_snps_schMan_panel.vcf
        #442,247 snp loci across all (sman homologous) chromosomes

#get only (sman homolgous) autosomal snps
vcftools \
    --vcf cohort_snps_schMan_panel.vcf \
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
    >cohort_snps_schMan_autosomal_panel.vcf
    #339,022 snp loci across all (sman homologous) chromosomes

#filter autosomal snps in LD
plink \
    --vcf cohort_snps_schMan_autosomal_panel.vcf \
    --allow-extra-chr \
    --indep-pairwise 25 5 0.20 \
    --out cohort_snps_schMan_autosomal_panel_LD-25-5-2
    #Pruning complete.  238702 of 339022 variants removed (100,320 remaining).

vcftools \
    --vcf cohort_snps_schMan_autosomal_panel.vcf \
    --exclude cohort_snps_schMan_autosomal_panel_LD-25-5-2.prune.out \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_schMan_autosomal_panel_LD-25-5-2.vcf
    #After filtering, kept 105 out of 105 Individuals
    #After filtering, kept 100320 out of a possible 339022 Sites


