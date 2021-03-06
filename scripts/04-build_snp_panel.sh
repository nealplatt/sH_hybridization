#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# 04-build_snp_panel.sh - after genotyping...the vcf files are filtered in a way
#       that generates multiple output files, for example SNPs on autosomal
#       chromosomes, vs. SNPs on autosomal chromosomes with LD fitlering

# Uses a combination of conda environments and singularity

cd $RESULTS_DIR/build_snp_panel

#outgroup identification
# ID        Species
#--------------------
#ERR119622  bovis
#ERR103048  bovis
#ERR310937  curassoni
#ERR119623  curassoni
#ERR084970  haematobium
#ERR037800  haematobium
#SRR433865  haematobium
#ERR103051  bovis
#ERR119612  guineensis
#ERR119613  intercalatum
#ERR310940  margrebowiei
#ERR539850  guineensis
#ERR539851  mattheei
#ERR539852  guineensis
#ERR539853  bovis
#ERR539854  intercalatum 
#ERR539855  mattheei
#ERR539856  intercalatum 
#ERR539857  mattheei

ln -s $SNP_DIR/cohort_filtered_bqsr-1.vcf
ln -s $SNP_DIR/cohort_filtered_bqsr-1.vcf.idx

${ENVIRONMENTS["SINGULARITY"]} \
    gatk SelectVariants \
        -V cohort_filtered_bqsr-1.vcf \
        -select-type SNP \
        -O cohort_snps.vcf \
        -R $HAE_GENOME
    #** 617,668 SNPs **#


#annotate snps names
bcftools annotate \
    --set-id +'%CHROM\:%POS' \
    cohort_snps.vcf \
    >cohort_snps_annotated.vcf

#filter loci with 85% of samples genotyped
vcftools \
    --vcf cohort_snps_annotated.vcf \
    --missing-site \
     --stdout \
    | sed 1d \
    | awk '{if ($6>0.15) print $1"\t"$2}' \
    >remove_site.list
    #removing 33,255 sites

vcftools \
    --vcf cohort_snps_annotated.vcf \
    --exclude-positions remove_site.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_85perSite.vcf
    #584,413 SNPs
#-------------------------

#filter individuals not genotyped at 20% of sites.
#   80% genotyping rate
vcftools \
    --vcf cohort_snps_85perSite.vcf \
    --missing-indv \
    --stdout \
    >indiv_missing_table.tsv

cat indiv_missing_table.tsv \
    | awk '$5 >0.20 {{print $1}}' \
    | sed 1d \
    >remove_indiv.list
    #ERR084970 S. haem (egypt) added manually; three samples from same seq library
    #SRR433865 S. haem (egypt) added manually; three samples from same seq library
    #ERR119622 S. bovis
    #ERR119623 S. curassoni
    #ERR539850 S. guineensis
    #ERR539851 S. bovis
    #ERR539852 S. guineensis
    #ERR539853 S. bovis
    #ERR539854 S. intercalatum
    #ERR539856 S. intercalatum
    #Sh.TZ_PEM0075.1
    #Sm.BR_0447.1
    #Sm.BR_1278.1
    #Sm.BR_2039.1
    #Sh.NE_Doki_029.1
    #Sh.NE_Kar_241.2

vcftools \
    --vcf cohort_snps_85perSite.vcf \
    --remove remove_indiv.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_80perIndiv.vcf
    #After filtering, kept 102 out of 118 Individuals
    #After filtering, kept 58,4413 out of a possible 58,4413 Sites


#-------------------------
#only keep bi-allelic sites
vcftools \
    --vcf cohort_snps_80perIndiv.vcf \
    --min-alleles 2 \
    --max-alleles 2 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_biallelic.vcf
    #After filtering, kept 102 out of 102 Individuals
    #After filtering, kept 549,249 out of a possible 584,413 Sites


#now get snps that are single locus in mansoni genome.
#lift over to smansoni coordinates

#convert to bed
vcf2bed \
    --do-not-sort \
    --max-mem=2G \
    <cohort_snps_biallelic.vcf \
    | awk '{print $1"\t"$2"\t"$3"\t"$4}' \
    >cohort_snps_biallelic.bed
    #549,249 remaining snps

#lift to sman coordinates
${ENVIRONMENTS["SINGULARITY"]} \
    /usr/software/progressiveCactus/submodules/hal/bin/halLiftover  \
        $RESULTS_DIR/wga/schMan7_vs_schMan1.hal \
        schHae_v1 \
        cohort_snps_biallelic.bed \
        schMan_v7 \
        cohort_snps_schMan.bed
        #628,324 remaining snps (increase is due to multilocus snps in mansoni genome)

#now lift the vcf to the sman coords (this will remove snp loci not aligned to sman)
python $WORK_DIR/scripts/lift_over_vcf.py \
    cohort_snps_schMan.bed \
    cohort_snps_biallelic.vcf \
    cohort_snps_schMan.vcf
    #531,201 remaining snps

#find loci that map to only a single region
awk '{print $4}' cohort_snps_schMan.bed \
    | sort \
    | uniq -c \
    | awk '$1>1 {print $2}' \
    >cohort_snps_schMan_remove.list
    #48,776 multiloci snps to remove


#now look for multiple snps aligned to a single locus
cut -f1,2,3 cohort_snps_schMan.bed \
    | sort \
    | uniq -c \
    | awk '$1>1 {print $2"\t"$3}' \
    >cohort_snps_schMan_remove.pos
    #4,823 loci annotated wtih multiple snps

#and put them in a vcf
vcftools \
    --vcf cohort_snps_schMan.vcf \
    --exclude-positions cohort_snps_schMan_remove.pos \
    --exclude cohort_snps_schMan_remove.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_schMan_singleLocus.vcf
    #After filtering, kept 102 out of 102 Individuals
    #After filtering, kept 482,406 out of a possible 531,201 Sites


#final check in vcf for multi-position snps
#check by snp_name
grep -v "#" cohort_snps_schMan_singleLocus.vcf \
    | cut -f3 \
    | sort \
    | uniq -c \
    | awk '{if ($1 > 1) print $2"\t"$3}' \
    >remove.list
    #0 redundant loci  

#check by pos
grep -v "#" cohort_snps_schMan_singleLocus.vcf \
    | cut -f1,2 \
    | sort \
    | uniq -c \
    | awk '{if ($1 > 1) print $2"\t"$3}' \
    >remove.pos
    #101 redundant loci

vcftools \
    --vcf cohort_snps_schMan_singleLocus.vcf \
    --exclude-positions remove.pos \
    --recode \
    --recode-INFO-all \
    --stdout \
    >tmp.vcf

mv tmp.vcf cohort_snps_schMan_singleLocus.vcf
#After filtering, kept 102 out of 102 Individuals
#After filtering, kept 482,204 out of a possible 482,406 Sites

#**********************************************************************************
#fix header and sort
grep -v "contig=<ID=" cohort_snps_schMan_singleLocus.vcf \
    >headerless.vcf

#add contigs for sman to header
${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk SelectVariants \
        -R $MAN_GENOME \
        -V headerless.vcf \
        -O header.vcf
        #482,204 loci remaining

#may need to be run on high mem nodes
${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk SortVcf \
        -R $MAN_GENOME \
        -I header.vcf \
        -O cohort_snps_schMan_panel.vcf
        ##482,204 snp loci across all (sman homologous) chromosomes

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
#After filtering, kept 102 out of 102 Individuals
#After filtering, kept 370770 out of a possible 482204 Sites

#filter all snps for maf
vcftools \
    --vcf cohort_snps_schMan_autosomal_panel.vcf \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_schMan_autosomal_maf05_panel.vcf
#After filtering, kept 102 out of 102 Individuals
#After filtering, kept 69,957 out of a possible 370770 Sites


#filter autosomal snps in LD
plink \
    --vcf cohort_snps_schMan_autosomal_maf05_panel.vcf\
    --allow-extra-chr \
    --indep-pairwise 25 5 0.20 \
    --out cohort_snps_schMan_autosomal_maf0_panel_LD
#Pruning complete.  64075 of 69957 variants removed.



vcftools \
    --vcf cohort_snps_schMan_autosomal_maf05_panel.vcf \
    --exclude cohort_snps_schMan_autosomal_maf0_panel_LD.prune.out \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_schMan_autosomal_maf05_LD.vcf
    #After filtering, kept 102 out of 102 Individuals
    #After filtering, kept 5,882 out of a possible 64,075 Sites


#once all filtering is complete there are three major autosomal files
#-------------------------------------------------------------------------------
#cohort_snps_schMan_autosomal_panel.vcf         all autosomal snps
#cohort_snps_schMan_autosomal_maf05_panel.vcf   all autosomal snps, maf filtered
#cohort_snps_schMan_autosomal_maf05_LD.vcf      all autosomal snps, maf filtered, ld filterd

#for ease...creating soft links
ln -s cohort_snps_schMan_autosomal_panel.vcf auto.vcf
ln -s cohort_snps_schMan_autosomal_maf05_panel.vcf auto_maf.vcf
ln -s cohort_snps_schMan_autosomal_maf05_LD.vcf auto_maf_ld.vcf

#and sample lists for downstream
cat build_snp_panel/auto_maf.vcf \
    | grep "#" \
    | tail -n 1 \
    | cut -f10- \
    | sed 's/\t/\n/gi' \
    | grep Sh.NE \
    >niger.list 


cat build_snp_panel/auto_maf.vcf \
    | grep "#" \
    | tail -n 1 \
    | cut -f10- \
    | sed 's/\t/\n/gi' \
    | grep Sh.TZ \
    >tz.list 



#presumed sex chromosome snps
vcftools \
    --vcf cohort_snps_schMan_panel.vcf \
    --chr SM_V7_ZW \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_schMan_ZW_panel.vcf
    #After filtering, kept 111,343 out of a possible 482,204 Sites

vcftools \
    --vcf cohort_snps_schMan_ZW_panel.vcf \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_schMan_ZW_maf05_panel.vcf
    #After filtering, kept 16940 out of a possible 111,343 Sites

#cohort_snps_schMan_ZW_panel.vcf                all ZW snps
#cohort_snps_schMan_ZW_maf05_panel.vcf          all ZW snps, maf filtered
ln -s cohort_snps_schMan_ZW_maf05_panel.vcf zw_maf.vcf
ln -s cohort_snps_schMan_ZW_panel.vcf zw.vcf


#haem
#haem - maf
#haem - LD
vcftools \
    --vcf cohort_snps_biallelic.vcf \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_schHae_maf05_panel.vcf
#After filtering, kept 102 out of 102 Individuals
#After filtering, kept 102020 out of a possible 549249 Sites


#filter autosomal snps in LD
plink \
    --vcf cohort_snps_schHae_maf05_panel.vcf\
    --allow-extra-chr \
    --indep-pairwise 25 5 0.20 \
    --out cohort_snps_schHae_maf05_panel_LD
    #Pruning complete.  90302 of 102020 variants removed.

vcftools \
    --vcf cohort_snps_schHae_maf05_panel.vcf \
    --exclude cohort_snps_schHae_maf05_panel_LD.prune.out \
    --recode \
    --recode-INFO-all \
    --stdout \
    >cohort_snps_schHae_maf05_LD.vcf
    #After filtering, kept 102 out of 102 Individuals
    #After filtering, kept 11,718 out of a possible 102,020 Sites

#cohort_snps_biallelic.vcf                  all biallelic alleles, haem coords
#cohort_snps_schHae_maf05_panel.vcf         bi-alleles, maf filtered, haem coords
#cohort_snps_schHae_maf05_LD.vcf     bi-alleles, maf and LD filtered, haem coords

ln -s cohort_snps_biallelic.vcf haem.vcf
ln -s cohort_snps_schHae_maf05_panel.vcf haem_maf.vcf
ln -s cohort_snps_schHae_maf05_LD.vcf haem_maf_ld.vcf


