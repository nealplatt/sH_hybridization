source activate snp_calling
source /master/nplatt/schisto_hybridization/scripts/set_env.sh

cd $RESULTS_DIR

mkdir whats_in_the_missing_regions

cd whats_in_the_missing_regions

cp ../niger.list .
cp ../tz.list .

#get the vcf for the missing region

vcftools \
    --vcf ../build_snp_panel/auto_maf.vcf \
    --keep niger.list \
    --keep tz.list \
    --chr SM_V7_4 \
    --from-bp 18000000 \
    --to-bp 20000000 \
    --recode \
    --stdout \
    >missing_region.vcf 

#scrolling through the data to find big breaks
grep -v "#" missing_region.vcf | cut -f1-3 | less

#...
#...
#SM_V7_4 18886348        KL252404.1:13951
#SM_V7_4 19828369        KL251477.1:20130
#...
#...

#break from SM_V7_4:18,886,348-19,828,369 ... ~1Mbp.

#are there genes in this region in Smansoni genome
#make bed file
echo "SM_V7_4 18886348  19828369" >missing_region.bed

#find genes that intersect this region
bedtools intersect -wb -a missing_region.bed -b ../../data/Sm_v7.0.bed | cut -f7 | cut -f1 -d"." | sort | uniq
#Smp_127050
#Smp_137250
#Smp_168430
#Smp_204920
#Smp_205780
#Smp_316180
#Smp_316200
#Smp_316210
#Smp_316220
#Smp_316230
#Smp_316240
#Smp_316250
#Smp_336130
#Smp_337070
#Smp_340550
#Smp_345820
#
# 17 of them

# so one of three things happened
#1) the wga didn't lift to these regions (or was filtered)
#2) we didn't design probes for these genes
#3) this region is missing in the haem (or bovis?) genome

# DID WE LIFT TO THIS REGION?
${ENVIRONMENTS["TITAN SINGULARITY"]} \
    /usr/software/progressiveCactus/submodules/hal/bin/halLiftover  \
        $RESULTS_DIR/wga/schMan7_vs_schMan1.hal \
        schMan_v7 \
        missing_region.bed \
        schHae_v1 \
        missing_region_schHae.bed

#looks like a hodge podge
#...
#KL251477.1      2983    3103
#KL251477.1      661     887
#KL252504.1      1298    1560
#AMPZ01027874.1  1531    1793
#KL252504.1      1275    1298
#AMPZ01027874.1  1793    1816
#KL251477.1      642     661
#KL251477.1      640     642
#KL251477.1      539     640
#KL252504.1      1240    1275
#AMPZ01027874.1  1816    1851
#...

#did we map any snps to these regions
bedtools intersect \
    -wo \
    -a missing_region_schHae.bed \
    -b ../build_snp_panel/cohort_snps_biallelic.bed \
    >snps_in_missing_region.bed

cut -f7 snps_in_missing_region.bed | sort | uniq | wc -l
#247 snps in this region

#where were they lost
head snps_in_missing_region.bed 
#...
#KL252404.1      13951   14260   KL252404.1      13981   13982   KL252404.1:13982        1
#KL251913.1      13293   13556   KL251913.1      13324   13325   KL251913.1:13325        1
#KL251913.1      13293   13556   KL251913.1      13341   13342   KL251913.1:13342        1
#KL251913.1      13293   13556   KL251913.1      13342   13343   KL251913.1:13343        1
#KL251913.1      13293   13556   KL251913.1      13352   13353   KL251913.1:13353        1
#KL251913.1      13293   13556   KL251913.1      13353   13354   KL251913.1:13354        1
#KL251913.1      13293   13556   KL251913.1      13361   13362   KL251913.1:13362        1
#KL251913.1      13293   13556   KL251913.1      13374   13375   KL251913.1:13375        1
#KL251913.1      13293   13556   KL251913.1      13377   13378   KL251913.1:13378        1
#KL251913.1      13293   13556   KL251913.1      13380   13381   KL251913.1:13381        1
#...

#find loci that map to only a single region
awk '{print $7}' snps_in_missing_region.bed \
    | sort \
    | uniq -c \
    | awk '{if ($1==1) print $2}' \
    >single_locus_snps.list

#only 8 of the 247 snps that lifted to this region are single locus snps
#KL251477.1:17244
#KL251477.1:17259
#KL251477.1:17274
#KL251477.1:17304
#KL251477.1:17328
#KL251477.1:17333
#KL251477.1:20130
#KL252404.1:13982

#so what happened to these 8 (6 of them are within 100bp of one another)
grep -f single_locus_snps.list snps_in_missing_region.bed >single_locus_snps.bed
#KL252404.1      13951   14260   KL252404.1      13981   13982   KL252404.1:13982        1
#KL251477.1      17191   17367   KL251477.1      17243   17244   KL251477.1:17244        1
#KL251477.1      17191   17367   KL251477.1      17258   17259   KL251477.1:17259        1
#KL251477.1      17191   17367   KL251477.1      17273   17274   KL251477.1:17274        1
#KL251477.1      17191   17367   KL251477.1      17303   17304   KL251477.1:17304        1
#KL251477.1      17191   17367   KL251477.1      17327   17328   KL251477.1:17328        1
#KL251477.1      17191   17367   KL251477.1      17332   17333   KL251477.1:17333        1
#KL251477.1      20048   20130   KL251477.1      20129   20130   KL251477.1:20130        1
# Are these positions with only 1 snp mapped to them?
cut -f4,5,6 single_locus_snps.bed >single_locus_snps_pos.list



cut -f4,5,6 snps_in_missing_region.bed \
    | sort \
    | uniq -c \
    | awk '$1>1 {print $2"\t"$3}' \
    >multiple_snps_to_same_pos.list

grep -f single_locus_snps_pos.list multiple_snps_to_same_pos.list
#nada

single_locus_snps.list
#where they removed due to low genotyping rate?
#generate list of snps removed due to low gt
vcftools \
    --vcf ../build_snp_panel/cohort_snps_annotated.vcf \
    --snps single_locus_snps.list \
    --missing-site \
    --stdout
#CHR             POS     N_DATA  N_GENOTYPE_FILTERED     N_MISS  F_MISS
#KL251477.1      17244   236     0                       20      0.0847458
#KL251477.1      17259   236     0                       18      0.0762712
#KL251477.1      17274   236     0                       18      0.0762712
#KL251477.1      17304   236     0                       18      0.0762712
#KL251477.1      17328   236     0                       20      0.0847458
#KL251477.1      17333   236     0                       14      0.059322
#KL251477.1      20130   236     0                       14      0.059322
#KL252404.1      13982   236     0                       32      0.135593

#what about maf?
vcftools \
    --vcf ../build_snp_panel/cohort_snps_annotated.vcf \
    --snps single_locus_snps.list \
    --freq \
    --stdout
#CHROM           POS     N_ALLELES   N_CHR   MAJ:FREQ        MIN:FREQ
#KL251477.1      17244   2           216     G:0.981481      A:0.0185185
#KL251477.1      17259   2           218     G:0.981651      A:0.0183486
#KL251477.1      17274   2           218     A:0.981651      T:0.0183486
#KL251477.1      17304   2           218     G:0.981651      A:0.0183486
#KL251477.1      17328   2           216     C:0.981481      T:0.0185185
#KL251477.1      17333   2           222     G:0.995495      A:0.0045045
#KL251477.1      20130   2           222     A:0.671171      G:0.328829
#KL252404.1      13982   2           204     G:0.995098      A:0.00490196

#this would filter out all except:
#KL251477.1:20130

#where was it removed
grep KL251477.1:20130 ../build_snp_panel/cohort_snps_schMan_singleLocus.vcf
#its here


grep KL251477.1:20130 ../build_snp_panel/cohort_snps_schMan_panel.vcf
SM_V7_4 19828369
#its here as well

#im an idiot - the span of the gap is (as noted above) ..
#...
#...
#SM_V7_4 18886348        KL252404.1:13951
#SM_V7_4 19828369        KL251477.1:20130
#...
#...

#this snp is in the final dataset and marks the end of the gap.

#all snps in gap accounted for...
#239/247 are snps that map to more than one place in the schMan genome
#008/247 are snps with MAF < 0.05
#001/247 is the ending snp in the gap

