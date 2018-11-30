#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# filter_whole_genome_snps.sh - filter genotyped sites from vcf files of whole
#   genome resequencing data. 

# Uses a conda node with 125gb of mem

#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd /master/nplatt/schisto_hybridization/results/build_whole_genome_snp_panel

ln -s ../call_whole_genome_snps/shae-sbov_whole-genome_filtered_geno-qual.vcf
ln -s ../call_whole_genome_snps/shae-sbov_whole-genome_filtered_geno-qual.vcf.idx

singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    gatk SelectVariants \
        -V shae-sbov_whole-genome_filtered_geno-qual.vcf \
        -select-type SNP \
        -O wg_snps.vcf \
        -R $HAE_GENOME
        #9,507,098 SNPs    

#annotate snps names so that each is uniq
bcftools annotate \
    --set-id +'%CHROM\:%POS' \
    wg_snps.vcf \
    >wg_snps_annotated.vcf

#filter individuals not genotyped at 20% of sites (80% genotyping rate)
vcftools \
    --vcf wg_snps_annotated.vcf \
    --missing-indv \
    --stdout \
    >wg_indiv_missing_table.tsv


#INDV            N_DATA      N_MISS      F_MISS
#ERR103048       9507098     380717      0.0400456
#Sh_DaiCP_276_1  9507098     476137      0.0500823
#Sh_Dai_044_1    9507098     232801      0.0244871
#Sh_Kar_001_1    9507098     395280      0.0415774
#Sh_Kar_37_2     9507098     109742      0.0115432
#Sh_Lata_078_1   9507098     197995      0.020826
#Sh_PEM_0103_1   9507098     269283      0.0283244
#Sh_PEM_0104_1   9507098     830355      0.0873405
#Sh_PEM_0130_1   9507098     204652      0.0215262
#Sh_Tiag_272_1   9507098     177763      0.0186979
#Sh_UNG_0038_1   9507098     668614      0.0703279
#Sh_UNG_0099_1   9507098     197688      0.0207937
#Sh_UNG_0121_1   9507098     228131      0.0239959

#all samples genotyped at more than 80% of sites...none to exclude


#------------------------------------------------------------------------------

#filter loci with 80% of samples genotyped (11 of 13)
vcftools \
    --vcf wg_snps_annotated.vcf \
    --missing-site \
     --stdout \
    | sed 1d \
    | awk '{if ($6>0.20) print $1"\t"$2}' \
    >wg_remove_site.list
    #428,142 sites to remove

vcftools \
    --vcf wg_snps_annotated.vcf \
    --exclude-positions wg_remove_site.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >wg_snps_80perSite.vcf
    #After filtering, kept 13 out of 13 Individuals
    #After filtering, kept 9,078,956 out of a possible 9,507,098 Sites
#-------------------------


#-------------------------
#only keep bi-allelic sites
vcftools \
    --vcf wg_snps_80perSite.vcf \
    --min-alleles 2 \
    --max-alleles 2 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >wg_snps_biallelic.vcf
    #After filtering, kept 13 out of 13 Individuals
    #After filtering, kept 8,884,113 out of a possible 9,078,956 Sites



#now get snps that are single locus in mansoni genome.
#lift over to smansoni coordinates

#convert to bed
vcf2bed \
    --do-not-sort \
    --max-mem=2G \
    <wg_snps_biallelic.vcf \
    | awk '{print $1"\t"$2"\t"$3"\t"$4}' \
    >wg_snps_biallelic.bed
    #8,884,113

#lift to sman coordinates
singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    /usr/software/progressiveCactus/submodules/hal/bin/halLiftover  \
        $RESULTS_DIR/wga/schMan7_vs_schMan1.hal \
        schHae_v1 \
        wg_snps_biallelic.bed \
        schMan_v7 \
        wg_snps_schMan.bed
        #5,769,419

#now lift the vcf to the sman coords (this will remove snp loci not aligned to sman)
python $WORK_DIR/scripts/lift_over_vcf.py \
    wg_snps_schMan.bed \
    wg_snps_biallelic.vcf \
    wg_snps_schMan.vcf
    #5,217,740

#find loci that map to only a single region
awk '{print $4}' wg_snps_schMan.bed \
    | sort \
    | uniq -c \
    | awk '$1>1 {print $2}' \
    >wg_snps_schMan_remove.list
    #306,249 multiloci snps to remove

#now look for multiple snps aligned to a single locus
cut -f1,2,3 wg_snps_schMan.bed \
    | sort \
    | uniq -c \
    | awk '$1>1 {print $2"\t"$3}' \
    >wg_snps_schMan_remove.pos
    #24,955 loci annotated wtih multiple snps

#and put them in a vcf
vcftools \
    --vcf wg_snps_schMan.vcf \
    --exclude-positions wg_snps_schMan_remove.pos \
    --exclude wg_snps_schMan_remove.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >wg_snps_schMan_singleLocus.vcf
#After filtering, kept 13 out of 13 Individuals
#After filtering, kept 4,910,992 out of a possible 5,217,740 Sites


#final check in vcf for multi-position snps
#check by snp_name
grep -v "#" wg_snps_schMan_singleLocus.vcf \
    | cut -f3 \
    | sort \
    | uniq -c \
    | awk '{if ($1 > 1) print $2"\t"$3}' \
    >remove.list
    #0 redundant loci  

#check by pos
grep -v "#" wg_snps_schMan_singleLocus.vcf \
    | cut -f1,2 \
    | sort \
    | uniq -c \
    | awk '{if ($1 > 1) print $2"\t"$3}' \
    >remove.pos
    #4,089 redundant loci

vcftools \
    --vcf wg_snps_schMan_singleLocus.vcf \
    --exclude-positions remove.pos \
    --exclude remove.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >tmp.vcf
    #After filtering, kept 13 out of 13 Individuals
    #After filtering, kept 4,902,735 out of a possible 4,910,992 Sites

mv tmp.vcf wg_snps_schMan_singleLocus.vcf

#**********************************************************************************
#fix header and sort
grep -v "contig=<ID=" wg_snps_schMan_singleLocus.vcf \
    >headerless.vcf

#add contigs for sman to header
singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    gatk SelectVariants \
        -R $MAN_GENOME \
        -V headerless.vcf \
        -O header.vcf


#may need to be run on high mem nodes
singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
    gatk SortVcf \
        -R $MAN_GENOME \
        -I header.vcf \
        -O wg_snps_schMan_panel.vcf
        ##  snp loci across all (sman homologous) chromosomes

#get only (sman homolgous) autosomal snps
vcftools \
    --vcf wg_snps_schMan_panel.vcf \
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
    >wg_snps_schMan_autosomal_panel.vcf
#After filtering, kept 13 out of 13 Individuals
#After filtering, kept 3,804,711 out of a possible 4,902,735 Sites

#wg_SM_V7_1     1,190,560
#wg_SM_V7_2       611,799
#wg_SM_V7_3       579,691
#wg_SM_V7_4       657,189
#wg_SM_V7_5       229,993
#wg_SM_V7_6       335,066
#wg_SM_V7_7       200,413


#impute/phase with beagle
grep -v "#" wg_snps_schMan_autosomal_panel.vcf  \
    | awk '{printf "%s\t%s\t%.6f\t%s\n", $1, $3, $2/287000, $2}' \
    >auto.map

#create a list of sample to imput genotypes - in this analysis keeping all
#  haem samples and then the "good" bovis sample
grep "#" wg_snps_schMan_autosomal_panel.vcf \
    | tail -n 1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >samples.list

#impute on each seperatley
for i in $(seq 1 7); do

    CHR=SM_V7_"$i"

    #extract autosome specific vcf for the samples of interest
    vcftools \
        --vcf wg_snps_schMan_autosomal_panel.vcf \
        --chr $CHR \
        --keep samples.list \
        --recode \
        --stdout \
        >wg_$CHR".vcf"

    sed -i 's/,assembly=schMan_v7.fa//gi' wg_$CHR.vcf 

    grep $CHR auto.map >wg_$CHR".map" 

    beagle \
        gt=wg_$CHR".vcf" \
        out=wg_$CHR"_beagle" \
        map=wg_$CHR".map" \
        nthreads=4 \
        window=300 \
        overlap=30 \
        niterations=100 \
        >>wg_$CHR"_beagle.log" 2>&1 &

done

wait 

gunzip wg_SM_V7_?_beagle.vcf.gz


#create one file of phased snps (for all autosomes)
vcfcombine \
    wg_SM_V7_1_beagle.vcf \
    wg_SM_V7_2_beagle.vcf \
    wg_SM_V7_3_beagle.vcf \
    wg_SM_V7_4_beagle.vcf \
    wg_SM_V7_5_beagle.vcf \
    wg_SM_V7_6_beagle.vcf \
    wg_SM_V7_7_beagle.vcf \
    >wg_auto_beagle.vcf

#filter autosomal snps in LD
plink \
    --vcf wg_auto_beagle.vcf \
    --allow-extra-chr \
    --indep-pairwise 25 5 0.20 \
    --out wg_auto_beagle_LD \
    --double-id
    #3,234,476 of 3,804,711 variants removed.


vcftools \
    --vcf wg_auto_beagle.vcf \
    --exclude wg_auto_beagle_LD.prune.out \
    --recode \
    --recode-INFO-all \
    --stdout \
    >wg_auto_beagle_LD.vcf
    #After filtering, kept 13 out of 13 Individuals
    #After filtering, kept 570235 out of a possible 3804711 Sites


