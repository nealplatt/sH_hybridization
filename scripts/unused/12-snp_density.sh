#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# 12-snp_density.sh - calculate the number of snps in various windows (as 
#        defined by bedtools makewindows

# Uses a conda to manage the enviroment

#Set up the environmentsource /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

mkdir $RESULTS_DIR/snp_density
cd $RESULTS_DIR/snp_density

#make windows
bedtools makewindows \
    -g ../../data/genome/Smansoni_v7.fa.fai \
    -w 100000 \
    -s 100000 \
    | grep '\<SM_V7_.\>' \
    >sman_100k-100k.windows

#get sample files
cp ../niger.list .
cp ../tz.list .

#group
vcftools \
    --vcf ../build_snp_panel/auto_maf.vcf \
    --keep samples.list \
    --recode \
    --stdout \
    >all.vcf

#niger
vcftools \
    --vcf ../build_snp_panel/auto_maf.vcf \
    --keep niger.list \
    --recode \
    --stdout \
    >niger.vcf

#bovis
vcftools \
    --vcf ../build_snp_panel/auto_maf.vcf \
    --indv ERR103048 \
    --recode \
    --stdout \
    >bovis.vcf

#tz
vcftools \
    --vcf ../build_snp_panel/auto_maf.vcf \
    --keep tz.list \
    --recode \
    --stdout \
    >tz.vcf

for TYPE in tz bovis niger all; do
    bcftools view \
        -m2 \
        -M2 \
        -v snps \
        $TYPE.vcf \
        >$TYPE"_biallelic.vcf"

    vcf2bed \
        --do-not-sort \
        --max-mem=2G \
        < $TYPE"_biallelic.vcf" \
        | awk '{print $1"\t"$2"\t"$3"\t"$4}' \
        >$TYPE"_biallelic.bed"
done

#intersect with windows (to get snp/kb counts)
bedtools intersect -c -a sman_100k-100k.windows -b all_biallelic.bed  >sman_100k-100k.counts


#combine into a single file (with a header)


