#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# 16-ld.sh - compare ld between schHae_V1 vs. lifted schMan_v7 coordinates.  

# Uses a conda to manage the enviroment

#Set up the environment
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir ld

cd ld

#generate a file that give R2 and bp between markers for schMan
mkdir schMan_ld
for i in $(seq 1 7); do
    CHR=SM_V7_"$i"

        #get vcf chromosome
        vcftools \
            --vcf ../build_snp_panel/auto_maf.vcf \
            --keep haem.list \
            --maf 0.1 \
            --chr $CHR \
            --recode \
            --stdout \
            >$CHR.vcf

        #convert to ped
        plink \
            --threads 6 \
            --vcf $CHR.vcf \
            --out $CHR \
            --recode12 \
            --allow-extra-chr

        #calc R2 for all snps on the chr
        plink \
            --threads 6 \
            --r2 \
            --file $CHR \
            --out $CHR \
            --allow-extra-chr \
            --ld-window-r2 0 \
            --ld-window 75000 \
            --ld-window-kb 100000


        awk '{print $0"\t"$5-$2}' $CHR.ld >schMan_ld/$CHR.ld
    rm $CHR.ped $CHR.map $CHR.nosex $CHR.log $CHR.vcf $CHR.ld
done



###############
# to do something similar with schHaem coordinates need to re-build vcf


#build new vcf file with schHaem coords
grep -v "#" ../build_snp_panel/auto_maf.vcf \
    | awk  '{print $3"\t"$0}' \
    | sed 's/:/\t/' \
    | awk  '{$3=$3":"$4;; print $0}' \
    | sed 's/ /\t/g' \
    | cut --complement -f4,5 \
    >vcf.entries

grep "#" ../build_snp_panel/auto_maf.vcf | grep -v "contig=<ID=" >vcf.sm_header   

cat vcf.sm_header vcf.entries >broke_header.vcf

#add contigs for sman to header
${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk SelectVariants \
        -R $HAE_GENOME \
        -V broke_header.vcf \
        -O header.vcf
        #75,704 loci remaining

#now sort
${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk SortVcf \
        -R $HAE_GENOME \
        -I header.vcf \
        -O auto_maf_schHae.vcf
        #75,704 snp loci across all (sman homologous) chromosomes

#make sure they are bi-allelic
bcftools view \
    -m2 \
    -M2 \
    -v snps \
    auto_maf_schHae.vcf \
    >auto_maf_schHae_bi.vcf

#remove low frequency snps
vcftools \
    --vcf auto_maf_schHae.vcf \
    --keep haem.list \
    --maf 0.1 \
    --recode \
    --stdout \
    >auto_maf_schHae_bi_haem.vcf

#create a list of contigs (have to do on a contig by contig basis (like chrs in 
# schMan)).
grep -v "#" auto_maf_schHae_bi_haem.vcf \
    | cut -f1 \
    | sort \
    | uniq -c \
    | awk '{ if ($1>1) print $2}' \
    >schHae_contigs.list                                                   

#generate a file that give R2 and bp between markers for schHae
mkdir schHae_ld
#calc ld for all schHae contigs
for CHR in $(cat schHae_contigs.list); do

    #get vcf chromosome
    vcftools \
        --vcf auto_maf_schHae.vcf \
        --keep haem.list \
        --maf 0.1 \
        --chr $CHR \
        --recode \
        --stdout \
        >$CHR".vcf"

    #convert to ped
    plink \
        --threads 6 \
        --vcf $CHR".vcf" \
        --out $CHR \
        --recode12 \
        --allow-extra-chr

    #calc R2 for all snps on the chr
    plink \
        --threads 6 \
        --r2 \
        --file $CHR \
        --out $CHR \
        --allow-extra-chr \
        --ld-window-r2 0 \
        --ld-window 75000 \
        --ld-window-kb 100000

    #calc distance
    awk '{print $0"\t"$5-$2}' $CHR.ld >schHae_ld/$CHR.ld
    #rm $CHR.ped $CHR.map $CHR.nosex $CHR.log $CHR.vcf $CHR.ld

done

#when finished, create csv files with r2 and physical distance
cat schHae_ld/*.ld | awk '{print $7","$8}' | grep -v "R2" >schHae_ld.csv
grep -v R2 schMan_ld/SM_V7_?.ld | awk '{print $7","$8 }' >schMan_ld.csv


#cleanup
tar -czf schHae_ld.tgz schHae_ld/
tar -czf schMan_ld.tgz schMan_ld/
rm -r schHae_ld/ schMan_ld
rm KL* AMP*

#plot schHae_ld.csv and schMan_ld.csv in R
