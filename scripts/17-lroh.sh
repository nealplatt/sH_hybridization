#LROH


# ABANDONED ANALYSIS DUE TO BIASED DISTRIBUTION OF EXOME/CAPTURE BAITS



#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir lroh

cd lroh

cp ../niger.list .
cp ../tz.list .

#identify regions of homozygosity in NE and TZ (and find regions that differ)
for i in $(seq 1 7); do
    CHR=SM_V7_"$i"
    for POP in niger tz; do
    
        vcftools \
            --vcf ../beagle/auto_beagle.vcf \
            --chr $CHR \
            --LROH \
            --keep $POP.list \
            --stdout \
            >$CHR"_"$POP.lroh

        sed '1d' $CHR"_"$POP.lroh >$CHR"_"$POP.bed
    done
done

cat SM_V7_?_tz.bed >tz.bed
cat SM_V7_?_niger.bed >niger.bed

bedtools genomecov \
    -d \
    -i tz.bed \
    -g ../../data/genome/schMan_v7.fa.fai \
    >tz_lroh_d.bed

bedtools genomecov \
    -d \
    -i niger.bed \
    -g ../../data/genome/schMan_v7.fa.fai \
    >niger_lroh_d.bed
      
#have homozygosity by site...need to plot

#shuffle and find expectation based on stdev

#how many homo bases
cat tz_lroh_d.bed | awk '{if ($3>0) print $3}' | wc -l
#254,497,324
#409,579,008
cat tz_lroh_d.bed \
    | awk '{if ($3>0) print $3}' \
    | sort -n \
    | uniq -c \
    | awk '{sum+=$1*$2} END {print sum}'
    
#this gives number of homozygous sites for the TZ population.
#to get the prob that a site will have 2 homozygous samples          

num homozygous sites        
