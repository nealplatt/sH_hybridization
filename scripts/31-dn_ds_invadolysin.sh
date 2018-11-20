source activate snp_calling
source /master/nplatt/schisto_hybridization/scripts/set_env.sh

cd $RESULTS_DIR

mkdir dnds

cd dnds

#mask genome for unsampled regions

#first convert all bams to bed and merge
for BAM in $(ls ../map_reads/exome/*_processed.bam); do
    
    SAMPLE=$(basename $BAM .bam)

    genomeCoverageBed \
        -bga \
        -split \
        -ibam $BAM \
        | bedtools \
            merge \
            -i - \
            >$SAMPLE.bed
done &


for BAM in $(ls ../map_reads/exome/*_processed.bam); do
    
    echo -e "    $BAM \\"
done &


cut -f1,2 ../../data/genome/schHae_v1.fa.fai >schHae_v1.fa.genome
awk '{print $1"\t0\t"$2"\t"$1}' schHae_v1.fa.genome >schHae_v1.fa.genome.bed

bedtools \
    multicov \
    -bed schHae_v1.fa.genome.bed \
    -bams \
    ../map_reads/exome/ERR037800_processed.bam \
    ../map_reads/exome/ERR084970_processed.bam \
    ../map_reads/exome/ERR103048_processed.bam \
    ../map_reads/exome/ERR103051_processed.bam \
    ../map_reads/exome/ERR119612_processed.bam \
    ../map_reads/exome/ERR119613_processed.bam \
    ../map_reads/exome/ERR119622_processed.bam \
    ../map_reads/exome/ERR119623_processed.bam \
    ../map_reads/exome/ERR310937_processed.bam \
    ../map_reads/exome/ERR310940_processed.bam \
    ../map_reads/exome/ERR539850_processed.bam \
    ../map_reads/exome/ERR539851_processed.bam \
    ../map_reads/exome/ERR539852_processed.bam \
    ../map_reads/exome/ERR539853_processed.bam \
    ../map_reads/exome/ERR539854_processed.bam \
    ../map_reads/exome/ERR539855_processed.bam \
    ../map_reads/exome/ERR539856_processed.bam \
    ../map_reads/exome/ERR539857_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-002.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-010.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-013.3_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-031.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-033.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-044.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-045.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-051.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-074.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-146.1_processed.bam \
    ../map_reads/exome/Sh.NE_DaiCP-233.1_processed.bam \
    ../map_reads/exome/Sh.NE_DaiCP-276.1_processed.bam \
    ../map_reads/exome/Sh.NE_DaiCP-277.2_processed.bam \
    ../map_reads/exome/Sh.NE_Doki-029.1_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-001.1_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-002.1_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-076.1_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-096.2_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-241.1_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-241.2_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-281.1_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-37.2_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-007.3_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-033.1_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-078.1_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-253.1_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-275.2_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-293.1_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-294.1_processed.bam \
    ../map_reads/exome/Sh.NE_LibTB-009.2_processed.bam \
    ../map_reads/exome/Sh.NE_LibTB-010.1_processed.bam \
    ../map_reads/exome/Sh.NE_LibTB-022.1_processed.bam \
    ../map_reads/exome/Sh.NE_LibTB-028.1_processed.bam \
    ../map_reads/exome/Sh.NE_LibTB-031.1_processed.bam \
    ../map_reads/exome/Sh.NE_NG-011.1_processed.bam \
    ../map_reads/exome/Sh.NE_NG-06.2_processed.bam \
    ../map_reads/exome/Sh.NE_NG-089.1_processed.bam \
    ../map_reads/exome/Sh.NE_NG-236.1_processed.bam \
    ../map_reads/exome/Sh.NE_Seb-076.1_processed.bam \
    ../map_reads/exome/Sh.NE_Seb-078.2_processed.bam \
    ../map_reads/exome/Sh.NE_Seb-081.2_processed.bam \
    ../map_reads/exome/Sh.NE_Tiag-272.1_processed.bam \
    ../map_reads/exome/Sh.NE_YK-029.2_processed.bam \
    ../map_reads/exome/Sh.NE_YK-069.1_processed.bam \
    ../map_reads/exome/Sh.NE_YK-099.2_processed.bam \
    ../map_reads/exome/Sh.NE_YK-248.2_processed.bam \
    ../map_reads/exome/Sh.NE_Youri-069.2_processed.bam \
    ../map_reads/exome/Sh.NE_Youri-091.3_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0063.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0075.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0076.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0079.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0089.2_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0094.2_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0099.2_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0103.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0104.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0106.2_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0108.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0110.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0114.3_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0115.4_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0120.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0125.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0126.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0127.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0128.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0130.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0133.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0139.2_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0145.3_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0154.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0157.3_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0166.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0171.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0006.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0038.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0076.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0077.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0078.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0087.2_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0089.3_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0092.3_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0099.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0102.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0111.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0117.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0121.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0125.3_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0127.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0129.2_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0134.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0137.3_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0139.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0142.2_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0146.1_processed.bam \
    ../map_reads/exome/Sm.BR_0447.1_processed.bam \
    ../map_reads/exome/Sm.BR_1278.1_processed.bam \
    ../map_reads/exome/Sm.BR_2039.1_processed.bam \
    ../map_reads/exome/SRR433865_processed.bam \
    >all_multicov.bed

#so i need to find all positions that have a bunch (n=?) of sh coverage
#  and atleast 3 of 6 for outgroup coverage.


#mask genome

#build cds for each gene
#make a list of each transcript
grep cds Sm_v7.0.bed | cut -f4 | sed 's/.cds//gi' | sort | uniq >cds.list
wc -l cds.list 
#14608 cds.list

mkdir cds_bed
grep cds Sm_v7.0.bed >all_cds.bed
#for each cds build a bed file that has the cds coordinates
while read TRANSCRIPT; do
    grep $TRANSCRIPT all_cds.bed \
        | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5*1000"\t"$6'} \
        >cds_bed/$TRANSCRIPT.bed
   
    #use this command to do a quick and dirty check that everything is in frame
    #awk '{SUM+=$3-$2} END {print SUM/3}' cds_bed/$TRANSCRIPT.bed

    #lift the coordinates to haematobium
    ${ENVIRONMENTS["TITAN SINGULARITY"]} \
        /usr/software/progressiveCactus/submodules/hal/bin/halLiftover  \
            $RESULTS_DIR/wga/schMan7_vs_schMan1.hal \
            schMan_v7 \
            cds_bed/$TRANSCRIPT.bed \
            schHae_v1 \
            cds_bed/$TRANSCRIPT"_schHae.bed"
 
    awk '{SUM+=$3-$2} END {print SUM/3}' cds_bed/$TRANSCRIPT"_schHae.bed"
    #check to see that things lifted over correctly
    #->haem_cleaned.bed

    #if so extract the fasta with haem bed file
    #->cds.fasta

    #depending on the orientation concatenate cds
    #->gene.fasta

done <head.list


#convert to necessary format for paml/codeml


################################################################################

#get only the invadolysin CDS
grep Smp_127030 all_cds.bed >Smp_127030_cds.bed

#convert from Sman to sHaem coords
${ENVIRONMENTS["TITAN SINGULARITY"]} \
    /usr/software/progressiveCactus/submodules/hal/bin/halLiftover  \
        $RESULTS_DIR/wga/schMan7_vs_schMan1.hal \
        schMan_v7 \
        Smp_127030_cds.bed \
        schHae_v1 \
        Smp_127030_cds_schHae.bed

#intersect the CDS with the probe coordinates
bedtools intersect -a ../../data/schHae_v1_probes.bed -b Smp_127030_cds_schHae.bed >Smp_127030_probed-cds_schHae.bed

#to calculate coverage across all probed regions...need to have single bp intervals
bedtools makewindows -b Smp_127030_probed-cds_schHae.bed  -w 1 -s 1 >Smp_127030_probed-cds_schHae_1bpWindows.bed

#calculate coverage across each base in probed regions of Smp_127030
bedtools \
    multicov \
    -bed Smp_127030_probed-cds_schHae_1bpWindows.bed \
    -bams \
    ../map_reads/exome/ERR037800_processed.bam \
    ../map_reads/exome/ERR103048_processed.bam \
    ../map_reads/exome/ERR103051_processed.bam \
    ../map_reads/exome/ERR119612_processed.bam \
    ../map_reads/exome/ERR119613_processed.bam \
    ../map_reads/exome/ERR310937_processed.bam \
    ../map_reads/exome/ERR310940_processed.bam \
    ../map_reads/exome/ERR539850_processed.bam \
    ../map_reads/exome/ERR539855_processed.bam \
    ../map_reads/exome/ERR539857_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-002.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-010.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-013.3_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-031.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-033.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-044.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-045.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-051.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-074.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-146.1_processed.bam \
    ../map_reads/exome/Sh.NE_DaiCP-233.1_processed.bam \
    ../map_reads/exome/Sh.NE_DaiCP-276.1_processed.bam \
    ../map_reads/exome/Sh.NE_DaiCP-277.2_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-001.1_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-002.1_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-076.1_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-096.2_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-241.2_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-281.1_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-37.2_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-007.3_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-033.1_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-078.1_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-253.1_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-275.2_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-293.1_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-294.1_processed.bam \
    ../map_reads/exome/Sh.NE_LibTB-009.2_processed.bam \
    ../map_reads/exome/Sh.NE_LibTB-010.1_processed.bam \
    ../map_reads/exome/Sh.NE_LibTB-022.1_processed.bam \
    ../map_reads/exome/Sh.NE_LibTB-028.1_processed.bam \
    ../map_reads/exome/Sh.NE_LibTB-031.1_processed.bam \
    ../map_reads/exome/Sh.NE_NG-011.1_processed.bam \
    ../map_reads/exome/Sh.NE_NG-06.2_processed.bam \
    ../map_reads/exome/Sh.NE_NG-089.1_processed.bam \
    ../map_reads/exome/Sh.NE_NG-236.1_processed.bam \
    ../map_reads/exome/Sh.NE_Seb-076.1_processed.bam \
    ../map_reads/exome/Sh.NE_Seb-078.2_processed.bam \
    ../map_reads/exome/Sh.NE_Seb-081.2_processed.bam \
    ../map_reads/exome/Sh.NE_Tiag-272.1_processed.bam \
    ../map_reads/exome/Sh.NE_YK-029.2_processed.bam \
    ../map_reads/exome/Sh.NE_YK-069.1_processed.bam \
    ../map_reads/exome/Sh.NE_YK-099.2_processed.bam \
    ../map_reads/exome/Sh.NE_YK-248.2_processed.bam \
    ../map_reads/exome/Sh.NE_Youri-069.2_processed.bam \
    ../map_reads/exome/Sh.NE_Youri-091.3_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0063.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0076.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0079.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0089.2_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0094.2_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0099.2_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0103.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0104.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0106.2_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0108.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0110.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0114.3_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0115.4_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0120.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0125.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0126.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0127.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0128.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0130.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0133.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0139.2_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0145.3_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0154.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0157.3_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0166.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0171.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0006.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0038.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0076.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0077.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0078.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0087.2_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0089.3_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0092.3_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0099.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0102.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0111.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0117.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0121.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0125.3_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0127.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0129.2_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0134.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0137.3_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0139.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0142.2_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0146.1_processed.bam \
    >Smp_127030_probed-cds_schHae_1bpWindows_multicov.tbl

#find probed regions of low coverage 
# to be included site must have:
#   85% Haem with >= 5 reads
#   gtd in bovis
#   gtd in curassoni
#   gtd in margrebowie
#   gtd in one other outgroup

#equal to lowest cov snp in invadolysin calls
vcftools \
    --vcf ../build_snp_panel/auto.vcf \
    --bed Smp_127030_cds.bed \
    --recode \
    --stdout >Smp_127030_cds.vcf

vcf2bed \
    --do-not-sort \
    --max-mem=2G \
    <Smp_127030_cds.vcf \
    | awk '{print $1"\t"$2"\t"$3"\t"$4}' \
    >Smp_127030_cds_snps.bed
    #60 remaining snps

${ENVIRONMENTS["TITAN SINGULARITY"]} \
    /usr/software/progressiveCactus/submodules/hal/bin/halLiftover  \
        $RESULTS_DIR/wga/schMan7_vs_schMan1.hal \
        schMan_v7 \
        Smp_127030_cds_snps.bed \
        schHae_v1 \
        Smp_127030_cds_snps_schHae.bed

bedtools \
    multicov \
    -bed Smp_127030_cds_snps_schHae.bed \
    -bams \
    ../map_reads/exome/ERR037800_processed.bam \
    ../map_reads/exome/ERR103048_processed.bam \
    ../map_reads/exome/ERR103051_processed.bam \
    ../map_reads/exome/ERR119612_processed.bam \
    ../map_reads/exome/ERR119613_processed.bam \
    ../map_reads/exome/ERR310937_processed.bam \
    ../map_reads/exome/ERR310940_processed.bam \
    ../map_reads/exome/ERR539850_processed.bam \
    ../map_reads/exome/ERR539855_processed.bam \
    ../map_reads/exome/ERR539857_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-002.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-010.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-013.3_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-031.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-033.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-044.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-045.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-051.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-074.1_processed.bam \
    ../map_reads/exome/Sh.NE_Dai-146.1_processed.bam \
    ../map_reads/exome/Sh.NE_DaiCP-233.1_processed.bam \
    ../map_reads/exome/Sh.NE_DaiCP-276.1_processed.bam \
    ../map_reads/exome/Sh.NE_DaiCP-277.2_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-001.1_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-002.1_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-076.1_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-096.2_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-241.2_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-281.1_processed.bam \
    ../map_reads/exome/Sh.NE_Kar-37.2_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-007.3_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-033.1_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-078.1_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-253.1_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-275.2_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-293.1_processed.bam \
    ../map_reads/exome/Sh.NE_Lata-294.1_processed.bam \
    ../map_reads/exome/Sh.NE_LibTB-009.2_processed.bam \
    ../map_reads/exome/Sh.NE_LibTB-010.1_processed.bam \
    ../map_reads/exome/Sh.NE_LibTB-022.1_processed.bam \
    ../map_reads/exome/Sh.NE_LibTB-028.1_processed.bam \
    ../map_reads/exome/Sh.NE_LibTB-031.1_processed.bam \
    ../map_reads/exome/Sh.NE_NG-011.1_processed.bam \
    ../map_reads/exome/Sh.NE_NG-06.2_processed.bam \
    ../map_reads/exome/Sh.NE_NG-089.1_processed.bam \
    ../map_reads/exome/Sh.NE_NG-236.1_processed.bam \
    ../map_reads/exome/Sh.NE_Seb-076.1_processed.bam \
    ../map_reads/exome/Sh.NE_Seb-078.2_processed.bam \
    ../map_reads/exome/Sh.NE_Seb-081.2_processed.bam \
    ../map_reads/exome/Sh.NE_Tiag-272.1_processed.bam \
    ../map_reads/exome/Sh.NE_YK-029.2_processed.bam \
    ../map_reads/exome/Sh.NE_YK-069.1_processed.bam \
    ../map_reads/exome/Sh.NE_YK-099.2_processed.bam \
    ../map_reads/exome/Sh.NE_YK-248.2_processed.bam \
    ../map_reads/exome/Sh.NE_Youri-069.2_processed.bam \
    ../map_reads/exome/Sh.NE_Youri-091.3_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0063.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0076.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0079.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0089.2_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0094.2_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0099.2_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0103.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0104.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0106.2_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0108.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0110.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0114.3_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0115.4_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0120.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0125.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0126.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0127.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0128.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0130.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0133.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0139.2_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0145.3_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0154.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0157.3_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0166.1_processed.bam \
    ../map_reads/exome/Sh.TZ_PEM0171.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0006.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0038.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0076.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0077.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0078.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0087.2_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0089.3_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0092.3_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0099.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0102.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0111.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0117.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0121.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0125.3_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0127.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0129.2_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0134.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0137.3_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0139.1_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0142.2_processed.bam \
    ../map_reads/exome/Sh.TZ_UNG0146.1_processed.bam \
    >Smp_127030_cds_snps_schHae_multicov.tbl

#so manually examining called snps vs all probed regions its apparent that we had
# good enough coverage to include all the probed regions in the dn ds caluculation
# roughly, the lowest called SNP had 11.5x coverage and all but 1 bp in the 
# entire probed regions of invadolysin had gt 11.5...and this region was like 10.5

#create bed intervals of unprobed Smp cds
bedtools intersect -v -a Smp_127030_cds_schHae.bed -b ../../data/schHae_v1_probes.bed  >Smp_127030_UNprobed-cds_schHae.bed

#need to build cds for invadolysin.
# first get invadolysin contig
samtools faidx ../../data/genome/schHae_v1.fa KL250964.1 >KL250964.1.fasta

#now mask unprobed regions
bedtools intersect -v -a Smp_127030_cds_schHae.bed -b ../../data/schHae_v1_probes.bed  >Smp_127030_UNprobed-cds_schHae.bed
bedtools maskfasta -fi KL250964.1.fasta -bed Smp_127030_UNprobed-cds_schHae.bed -fo KL250964.1_masked.fasta
samtools faidx KL250964.1_masked.fasta
${ENVIRONMENTS["TITAN SINGULARITY"]} gatk CreateSequenceDictionary -R KL250964.1_masked.fasta

#then fasta from reference of invadolysin snps
#...need to convert vcf from sman to shaem coords
schMan_to_schHae_vcf_coords.py Smp_127030_cds.vcf tmp.vcf 

Smp_127030_cds_schHae_coords.vcf
grep -v "contig=<ID=" tmp.vcf \
    >headerless.vcf

#add contigs for sman to header
${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk SelectVariants \
        -R $HAE_GENOME \
        -V headerless.vcf \
        -O header.vcf

#may need to be run on high mem nodes
${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk SortVcf \
        -R $HAE_GENOME \
        -I header.vcf \
        -O Smp_127030_cds_schHae_coords.vcf

#for each sample make a fasta
if [ -f Smp_127030_cds.fasta ]; then
    rm Smp_127030_cds.fasta
fi

for SAMPLE in $(cat ../samples.list); do
    
    #get the vcf
    vcftools \
        --vcf Smp_127030_cds_schHae_coords.vcf \
        --indv $SAMPLE \
        --recode \
        --recode-INFO-all \
        --stdout \
        >$SAMPLE"_Smp_127030_cds.vcf"

    #get the fasta
    gatk \
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

#now get the mansoni sequence and add to use as an outgroup
samtools faidx ../../data/genome/schMan_v7.fa SM_V7_4 >SM_V7_4.fasta

${ENVIRONMENTS["TITAN SINGULARITY"]} \
    /usr/software/progressiveCactus/submodules/hal/bin/halLiftover  \
        $RESULTS_DIR/wga/schMan7_vs_schMan1.hal \
        schHae_v1 \
        Smp_127030_UNprobed-cds_schHae.bed \
        schMan_v7 \
        Smp_127030_UNprobed-cds_schMan.bed

bedtools maskfasta -fi SM_V7_4.fasta -bed Smp_127030_UNprobed-cds_schMan.bed -fo SM_V7_4_masked.fasta

${ENVIRONMENTS["TITAN SINGULARITY"]} \
    /usr/software/progressiveCactus/submodules/hal/bin/halLiftover  \
        $RESULTS_DIR/wga/schMan7_vs_schMan1.hal \
        schHae_v1 \
        Smp_127030_cds_schHae.bed \
        schMan_v7 \
        Smp_127030_cds_schMan.bed

bedtools getfasta \
        -bed Smp_127030_cds_schMan.bed \
        -fi SM_V7_4_masked.fasta \
        -fo schMan_Smp_127030_cds_indiv.fasta

grep -v "^>" schMan_Smp_127030_cds_indiv.fasta \
        | awk 'BEGIN { ORS=""; print ">schMan\n" } { print }' > schMan_Smp_127030_combined_cds.fasta

    cat schMan_Smp_127030_combined_cds.fasta >>Smp_127030_cds.fasta
    echo >>Smp_127030_cds.fasta

#generate a tree with mansoni as the outgroup
raxmlHPC \
    -fa \
    -m GTRGAMMA \
    -p 12345 \
    -x 12345 \
    -#100 \
    -s Smp_127030_cds.fasta \
    -T 2 \
    -n 100boot \
    -o schMan

#add label to bovis + NE node manually in figtree

#convert all Ns to ?s in alignment and convert alignment to phylip 
# (two spaces between name and seq)

#remove " ' "s from the tree

#make seleciton and neutral control files
codeml neutral.ctl

#should i be doing this on haplotype data?
#run beagle on everything.
#calculate a recomb map.  1cm = 287Kb based on data from Criscione et al. (2009) avg
grep -v "#" ../build_snp_panel/auto.vcf  \
    | awk '{printf "%s\t%s\t%.6f\t%s\n", $1, $3, $2/287000, $2}' \
    >auto.map

#create a list of sample to imput genotypes - in this analysis keeping all
#  haem samples and then the "good" bovis sample
grep "#" ../build_snp_panel/auto.vcf \
    | tail -n 1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >samples.list

#impute on each seperatley
for i in $(seq 1 7); do

    CHR=SM_V7_"$i"

#    #extract autosome specific vcf for the samples of interest
#    vcftools \
#        --vcf ../build_snp_panel/auto.vcf \
#        --chr $CHR \
#        --recode \
#        --stdout >$CHR".vcf"

#    sed -i 's/,assembly=schMan_v7.fa//gi' $CHR.vcf

#    grep $CHR auto.map >$CHR".map"

    beagle \
        gt=$CHR".vcf" \
        out=$CHR"_beagle_all_samples" \
        map=$CHR".map" \
        nthreads=4 \
        window=300 \
        overlap=30 \
        niterations=250 >>$CHR.beaglelog 2>&1 &

done

wait 

gunzip SM_V7_*_beagle_all_samples.vcf.gz

#create one file of phased snps (for all autosomes)
vcfcombine \
    SM_V7_1_beagle_all_samples.vcf \
    SM_V7_2_beagle_all_samples.vcf \
    SM_V7_3_beagle_all_samples.vcf \
    SM_V7_4_beagle_all_samples.vcf \
    SM_V7_5_beagle_all_samples.vcf \
    SM_V7_6_beagle_all_samples.vcf \
    SM_V7_7_beagle_all_samples.vcf \
    >auto_beagle_all_samples.vcf
    #370,770

#now get the phased invadolysin haplotypes
 vcftools \
    --vcf auto_beagle_all_samples.vcf \
    --bed Smp_127030_cds.bed \
    --recode \
    --recode-INFO-all \
    --stdout \
    >Smp_127030_phased.vcf



#now get this vcf for phased invadolysin haplotypes
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

if [ -f Smp_127030_haplotypes_cds.fasta ]; then
    rm Smp_127030_haplotypes_cds.fasta
fi

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

#add schMan to the file
cat schMan_Smp_127030_combined_cds.fasta >>Smp_127030_haplotypes_cds.fasta
echo >>Smp_127030_haplotypes_cds.fasta

#get unique haplotypes
fasta_formatter -t -i Smp_127030_haplotypes_cds.fasta -o Smp_127030_haplotypes_cds.tab

sort -u -k2 Smp_127030_haplotypes_cds.tab | awk '{print ">"$1"\n"$2}' >Smp_127030_uniq_haplotypes_cds.fasta

#examined tree and the frequency of sequences to find a reduced subset for PAML 
# runs

grep -e A.ERR310940 \
    -e A.Sh.TZ_PEM0104.1 \
    -e A.ERR103051 \
    -e A.ERR103048 \
    -e A.Sh.NE_Tiag-272.1 \
    -e schMan \
    Smp_127030_haplotypes_cds.tab \
    | awk '{print ">"$1"\n"$2}' \
        >reduced_Smp_127030_uniq_haplotypes_cds.fasta

#find the best tree
raxmlHPC \
    -m GTRCAT \
    -p 12345 \
    -s reduced_Smp_127030_uniq_haplotypes_cds.fasta \
    -# 1000 \
    -n reduced_Smp_127030_uniq_haplotypes_cds \
    -o schMan

#clean up replicates
mkdir raxml_reduced_trees
mv RAxML_*.reduced_Smp_127030_uniq_haplotypes_cds.RUN.* raxml_reduced_trees/


#now start paml/codeml analyses by making selection and control files
# first get the best raxml tree and and indicate branch of interest
cat RAxML_bestTree.reduced_Smp_127030_uniq_haplotypes_cds
#((((A.ERR103051:0.00618220804948685539,A.Sh.TZ_PEM0104.1:0.00280672070595951260):0.00246171314196946042,A.ERR310940:0.00079230556832455656):0.00367147403398323602,(A.Sh.NE_Tiag-272.1:0.00060498047741618143,A.ERR103048:0.00057900965475889772):0.00217271227828123772):0.02950785716451633997,schMan:0.02950785716451633997);


#added branch label for PAML
echo "((((A.ERR103051:0.00618220804948685539,A.Sh.TZ_PEM0104.1:0.00280672070595951260):0.00246171314196946042,A.ERR310940:0.00079230556832455656):0.00367147403398323602,(A.Sh.NE_Tiag-272.1 #1:0.00060498047741618143,A.ERR103048:0.00057900965475889772)#1:0.00217271227828123772):0.02950785716451633997,schMan:0.02950785716451633997);" \
    >reduced_Smp_127030_uniq_haplotypes_cds.tree

# now convert fasta to phylip
fasta_formatter \
    -t \
    -i reduced_Smp_127030_uniq_haplotypes_cds.fasta \
    -o reduced_Smp_127030_uniq_haplotypes_cds.phy

#then add " 6 2262\n" to the beginning to make it a phylip formatted alignment

#convert all Ns to ?s in the sequence
sed -i 's/N/?/g' reduced_Smp_127030_uniq_haplotypes_cds.phy
#then fix a the ?E on the NE sample
sed -i 's/N?/NE/' reduced_Smp_127030_uniq_haplotypes_cds.phy
#fix tabs to spaces
sed -i 's/\t/    /' reduced_Smp_127030_uniq_haplotypes_cds.phy

#make the neutral and selection codeml control files and upload
codeml neutral.ctl
codeml selection.ctl


grep lnL reduced_Smp_127030_uh-branch_site_neutral.H0.mlc 
#lnL(ntime: 10  np: 14):  -2845.725768      +0.000000
grep lnL reduced_Smp_127030_uh-branch_site_selection.H1.mlc
#lnL(ntime: 10  np: 15):  -2845.103222      +0.000000

#LRT = 2*(-2845.103222 - -2845.725768) = 1.245092
chi2
#Chi-square critical values

#                                Significance level

# DF    0.9950   0.9750   0.9000   0.5000   0.1000   0.0500   0.0100   0.0010
#  1    0.0000   0.0010   0.0158   0.4549   2.7055   3.8415   6.6349  10.8276

#reject model with selection 
#now run with ALL uniq haplotypes-----------------------------------------------
#make the tree
raxmlHPC-PTHREADS \
        -T 2 \
        -m GTRCAT \
        -p 12345 \
        -# 1000 \
        -s Smp_127030_uniq_haplotypes_cds.fasta \
        -n Smp_127030_uniq_haplotypes_cds \
        -o schMan

#clean up a bit
mkdir raxml_all_uniq_haplotypes_trees
mv RAxML_*Smp_127030_uniq_haplotypes_cds* raxml_all_uniq_haplotypes_trees/

cat raxml_all_uniq_haplotypes_trees/RAxML_bestTree.Smp_127030_uniq_haplotypes_cds
#((A.ERR310940:0.00133451934108898816,((((A.ERR310937:0.00063928626593180151,A.Sh.NE_NG-011.1:0.00000100000050002909):0.00257661528568875432,A.ERR119612:0.00000100000050002909):0.00064065471620128247,((B.Sh.NE_Dai-013.3:0.00128854550202245474,(A.ERR103048:0.00064235093003460193,((A.Sh.NE_Kar-001.1:0.00064390292987796566,(A.ERR037800:0.00000100000050002909,A.Sh.NE_YK-029.2:0.00064491744825488509):0.00000100000050002909):0.00064390295564639786,A.Sh.NE_Lata-078.1:0.00000100000050002909):0.00000100000050002909):0.00064503434630353213):0.00064247399058681251,B.ERR119612:0.00000100000050002909):0.00000100000050002909):0.00390783835836038505,(A.Sh.TZ_PEM0103.1:0.00063638485263258643,(B.Sh.TZ_PEM0120.1:0.00063722371357710162,((((A.Sh.TZ_UNG0117.1:0.00000100000050002909,(A.Sh.TZ_UNG0139.1:0.00127094925457619102,(B.Sh.TZ_PEM0099.2:0.00000100000050002909,(A.Sh.TZ_PEM0133.1:0.00126873285883879398,((B.ERR103051:0.00126866320274240702,(A.ERR103051:0.00126896197669645774,A.ERR539855:0.00000100000050002909):0.00190872044923837223):0.00000100000050002909,(B.Sh.TZ_UNG0087.2:0.00000100000050002909,(B.Sh.TZ_PEM0103.1:0.00000100000050002909,(B.Sh.TZ_PEM0145.3:0.00064087518460384724,(((((A.Sh.TZ_PEM0099.2:0.00063772546284204606,A.Sh.TZ_UNG0134.1:0.00000100000050002909):0.00063606639400788645,A.Sh.TZ_PEM0139.2:0.00127830321049212498):0.00000100000050002909,(((B.Sh.TZ_UNG0129.2:0.00321759943031953719,(((B.Sh.TZ_PEM0171.1:0.00127245111283949314,A.Sh.TZ_PEM0063.1:0.00000100000050002909):0.00000100000050002909,A.Sh.TZ_UNG0111.1:0.00063465392391950383):0.00127404735379609566,B.Sh.TZ_PEM0126.1:0.00000100000050002909):0.00063514446530332769):0.00000100000050002909,A.Sh.TZ_PEM0114.3:0.00000100000050002909):0.00000100000050002909,B.Sh.TZ_PEM0133.1:0.00063519043260941411):0.00063603530497547388):0.00000100000050002909,A.Sh.TZ_PEM0154.1:0.00000100000050002909):0.00000100000050002909,A.Sh.TZ_UNG0125.3:0.00192060815928473847):0.00063695393764512983):0.00063467463078119637):0.00063447937791556610):0.00063365040618427440):0.00000100000050002909):0.00063540095162459905):0.00000100000050002909):0.00191099797591265039):0.00449467180776074327,(B.Sh.TZ_UNG0038.1:0.00063499030110272434,((A.Sh.TZ_PEM0079.1:0.00127841087972491228,(B.Sh.TZ_PEM0110.1:0.00516891239923836752,B.Sh.TZ_UNG0111.1:0.00000100000050002909):0.00063494227797514409):0.00000100000050002909,(B.Sh.TZ_PEM0079.1:0.00000100000050002909,(((((B.Sh.TZ_UNG0076.1:0.00063587060918523499,A.Sh.TZ_PEM0094.2:0.00000100000050002909):0.00063717516739843249,(B.Sh.TZ_UNG0117.1:0.00063469068797979773,A.Sh.TZ_PEM0076.1:0.00127420199109140836):0.00000100000050002909):0.00000100000050002909,B.Sh.TZ_PEM0114.3:0.00000100000050002909):0.00000100000050002909,B.Sh.TZ_PEM0139.2:0.00063552392336664316):0.00063733525493169072,A.Sh.TZ_UNG0087.2:0.00063578704845731245):0.00000100000050002909):0.00000100000050002909):0.00000100000050002909):0.00063577085734859925):0.00000100000050002909,B.Sh.TZ_UNG0146.1:0.00000100000050002909):0.00063663384739386401,A.Sh.TZ_PEM0145.3:0.00000100000050002909):0.00063906355174048407):0.00000100000050002909):0.00258673892333877772):0.00125809983233602274):0.03392093211114020901,schMan:0.03392093211114020901);

#add the foreground branch label
echo "((A.ERR310940:0.00133451934108898816,((((A.ERR310937:0.00063928626593180151,A.Sh.NE_NG-011.1:0.00000100000050002909):0.00257661528568875432,A.ERR119612:0.00000100000050002909):0.00064065471620128247,((B.Sh.NE_Dai-013.3:0.00128854550202245474,(A.ERR103048:0.00064235093003460193,((A.Sh.NE_Kar-001.1:0.00064390292987796566,(A.ERR037800:0.00000100000050002909,A.Sh.NE_YK-029.2:0.00064491744825488509):0.00000100000050002909):0.00064390295564639786,A.Sh.NE_Lata-078.1:0.00000100000050002909):0.00000100000050002909):0.00064503434630353213):0.00064247399058681251,B.ERR119612:0.00000100000050002909):0.00000100000050002909)#1:0.00390783835836038505,(A.Sh.TZ_PEM0103.1:0.00063638485263258643,(B.Sh.TZ_PEM0120.1:0.00063722371357710162,((((A.Sh.TZ_UNG0117.1:0.00000100000050002909,(A.Sh.TZ_UNG0139.1:0.00127094925457619102,(B.Sh.TZ_PEM0099.2:0.00000100000050002909,(A.Sh.TZ_PEM0133.1:0.00126873285883879398,((B.ERR103051:0.00126866320274240702,(A.ERR103051:0.00126896197669645774,A.ERR539855:0.00000100000050002909):0.00190872044923837223):0.00000100000050002909,(B.Sh.TZ_UNG0087.2:0.00000100000050002909,(B.Sh.TZ_PEM0103.1:0.00000100000050002909,(B.Sh.TZ_PEM0145.3:0.00064087518460384724,(((((A.Sh.TZ_PEM0099.2:0.00063772546284204606,A.Sh.TZ_UNG0134.1:0.00000100000050002909):0.00063606639400788645,A.Sh.TZ_PEM0139.2:0.00127830321049212498):0.00000100000050002909,(((B.Sh.TZ_UNG0129.2:0.00321759943031953719,(((B.Sh.TZ_PEM0171.1:0.00127245111283949314,A.Sh.TZ_PEM0063.1:0.00000100000050002909):0.00000100000050002909,A.Sh.TZ_UNG0111.1:0.00063465392391950383):0.00127404735379609566,B.Sh.TZ_PEM0126.1:0.00000100000050002909):0.00063514446530332769):0.00000100000050002909,A.Sh.TZ_PEM0114.3:0.00000100000050002909):0.00000100000050002909,B.Sh.TZ_PEM0133.1:0.00063519043260941411):0.00063603530497547388):0.00000100000050002909,A.Sh.TZ_PEM0154.1:0.00000100000050002909):0.00000100000050002909,A.Sh.TZ_UNG0125.3:0.00192060815928473847):0.00063695393764512983):0.00063467463078119637):0.00063447937791556610):0.00063365040618427440):0.00000100000050002909):0.00063540095162459905):0.00000100000050002909):0.00191099797591265039):0.00449467180776074327,(B.Sh.TZ_UNG0038.1:0.00063499030110272434,((A.Sh.TZ_PEM0079.1:0.00127841087972491228,(B.Sh.TZ_PEM0110.1:0.00516891239923836752,B.Sh.TZ_UNG0111.1:0.00000100000050002909):0.00063494227797514409):0.00000100000050002909,(B.Sh.TZ_PEM0079.1:0.00000100000050002909,(((((B.Sh.TZ_UNG0076.1:0.00063587060918523499,A.Sh.TZ_PEM0094.2:0.00000100000050002909):0.00063717516739843249,(B.Sh.TZ_UNG0117.1:0.00063469068797979773,A.Sh.TZ_PEM0076.1:0.00127420199109140836):0.00000100000050002909):0.00000100000050002909,B.Sh.TZ_PEM0114.3:0.00000100000050002909):0.00000100000050002909,B.Sh.TZ_PEM0139.2:0.00063552392336664316):0.00063733525493169072,A.Sh.TZ_UNG0087.2:0.00063578704845731245):0.00000100000050002909):0.00000100000050002909):0.00000100000050002909):0.00063577085734859925):0.00000100000050002909,B.Sh.TZ_UNG0146.1:0.00000100000050002909):0.00063663384739386401,A.Sh.TZ_PEM0145.3:0.00000100000050002909):0.00063906355174048407):0.00000100000050002909):0.00258673892333877772):0.00125809983233602274):0.03392093211114020901,schMan:0.03392093211114020901);" \
    >Smp_127030_uniq_haplotypes_cds.tree

#make the phylip file
fasta_formatter \
    -t \
    -i Smp_127030_uniq_haplotypes_cds.fasta \
    -o Smp_127030_uniq_haplotypes_cds.tab

cut -f1 Smp_127030_uniq_haplotypes_cds.tab >names
cut -f2 Smp_127030_uniq_haplotypes_cds.tab | sed 's/N/\?/g' >seqs

paste names seqs >Smp_127030_uniq_haplotypes_cds.phy

#then add " 50 2262\n" to the beginning to make it a phylip formatted alignment

#fix tabs to spaces
sed -i 's/\t/    /' Smp_127030_uniq_haplotypes_cds.phy

#make the ctl files
sed 's/reduced_//' \
    <neutral.ctl \
    >neutral_all_uniq_haps/neutral_all_uniq_haps.ctl
cp Smp_127030_uniq_haplotypes_cds.tree neutral_all_uniq_haps/
cp Smp_127030_uniq_haplotypes_cds.phy neutral_all_uniq_haps/

sed 's/reduced_//' \
    <selection.ctl \
    >selection_all_uniq_haps/selection_all_uniq_haps.ctl
cp Smp_127030_uniq_haplotypes_cds.tree selection_all_uniq_haps/
cp Smp_127030_uniq_haplotypes_cds.phy selection_all_uniq_haps/

sed 's/reduced_//' \
    <dnds.ctl \
    >dnds_all_uniq_haps/dnds_all_uniq_haps.ctl    
cp Smp_127030_uniq_haplotypes_cds.tree dnds_all_uniq_haps/
cp Smp_127030_uniq_haplotypes_cds.phy dnds_all_uniq_haps/

codeml neutral_all_uniq_haps.ctl
codeml selection_all_uniq_haps.ctl
codeml dnds_all_uniq_haps.ctl

#site-test
#neutral lnL(ntime: 98  np:102):  -3392.298001      +0.000000
#selection lnL(ntime: 98  np:103):  -3392.298006      +0.000000

#SM_V7_4 20033013        KL250964.1:40108 

