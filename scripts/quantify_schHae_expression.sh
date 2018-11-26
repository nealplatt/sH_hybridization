#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# quantify_schHae_expression.sh - try to estimate expression of invadolysin
#       gene(s) using RNA data from other experiments  

# Uses a conda and singularity to manage the enviroment; needs relativley high
#       mem compute (used Titan - 125Gb)

#Set up the environment
source activate rna_exp
source /master/nplatt/schisto_hybridization/scripts/set_env.sh


cd $RESULTS_DIR
mkdir haem_trans
cd haem_trans

mkdir logs

#===============================================================================
#Run	        Treatment	    InsertSize	    MBases	    Notes	           # 
#-------------------------------------------------------------------------------
#SRX3632881	    Eggs	        0	            2,159	    >1K eggs	    
#SRX3632879	    Adult_female	0	            1,846	    50 worms	    
#SRX3632877	    Adult_male	    0	            2,190	    50 worms	    

#----------------------------------------------------
#get the data from the SRA
for SRA_ACCESSION in SRX3632881 SRX3632879 SRX3632877; do
    fastq-dump --split-files --outdir $SRA_ACCESSION"_reads" --gzip $SRA_ACCESSION &
done


#was concerned about converting between assemblies...hit to many difficulties
# will lift at the end (even if manually)
####----------------------------------------------------
####to use the schMan gene names, need to lift coords between schman and schhae
###gff2bed <../../data/Sm_v7.0.gff >Sm_v7.0.bed

####creat a gff file that has been lifted between assemblies
###/opt/projects/singularity-2.4.2/bin/singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img \
###    /usr/software/progressiveCactus/submodules/hal/bin/halLiftover  \
###        --outBedVersion 12\
###        $RESULTS_DIR/wga/schMan7_vs_schMan1.hal \
###        schMan_v7 \
###        Sm_v7.0.bed \
###        schHae_v1 \
###        Sm_v7_0_to_schHae_v1.bed

####extract only the genes from the bed file
###grep Sm_v7_0_to_schHae_v1.bed


####----------------------------------------------------
#create a softlink to the haem genome
ln -s ../../data/genome/schHae_v1.fa

#build bowtie2 index
bowtie2-build  schHae_v1.fa schHae_v1

#map with tophat
for SRA_ACCESSION in SRX3632881 SRX3632879 SRX3632877; do

     tophat \
         --num-threads 12 \
         --output-dir $SRA_ACCESSION"_tophat" \
         schHae_v1 \
         $SRA_ACCESSION"_reads"/$SRA_ACCESSION"_1.fastq.gz" \
         $SRA_ACCESSION"_reads"/$SRA_ACCESSION"_2.fastq.gz"
done

#----------------------------------------------------
#cufflinks
# used to predict intron/exon structure and transcripts
for SAMPLE in SRX3632881 SRX3632879 SRX3632877; do

    cufflinks \
        --no-update-check \
        --num-threads 12 \
        --multi-read-correct \
        --output-dir $SAMPLE"_cufflinks" \
        --frag-bias-correct schHae_v1.fa \
        --library-type fr-unstranded \
        --upper-quartile-norm \
        $SAMPLE"_tophat"/accepted_hits.bam 

done

#----------------------------------------------------
#create a list of all the predicted transcripts
ls SRX*cufflinks/transcripts.gtf >assembly_GTF_list.txt

#then use cuffmerge to merge the predicted transcripts across all samples to
# generate a single/consensus/merged annotation
cuffmerge \
    --no-update-check \
    -o cuffmerge \
    --num-threads 12 \
    --ref-sequence schHae_v1.fa \
    assembly_GTF_list.txt 



#ran into errors switching to hisat2 - installed with conda
hisat2-build schHae_v1.fa schHae_v1

for SRA_ACCESSION in SRX3632881 SRX3632879 SRX3632877; do

    hisat2 \
        -q \
        --threads 12 \
        -x schHae_v1 \
        -1 $SRA_ACCESSION"_reads"/$SRA_ACCESSION"_1.fastq.gz" \
        -2 $SRA_ACCESSION"_reads"/$SRA_ACCESSION"_2.fastq.gz" \
        -S $SRA_ACCESSION"_tophat".sam

done


for SRA_ACCESSION in SRX3632881 SRX3632879 SRX3632877; do

    #sort sam file
    samtools view  -Sb $SRA_ACCESSION"_tophat.sam" \
        | samtools sort - $SRA_ACCESSION"_sorted"

    #use stringtie to assemble transcripts and get abundance
    stringtie \
        $SRA_ACCESSION"_sorted.bam" \
        -o $SRA_ACCESSION"_stringtie.gtf" \
        -A $SRA_ACCESSION"_stringtie_gene_abund.tab"

done

ls SRX*gtf >gtf.list

stringtie \
    --merge \
    -o schHaem_stringtie-merged.gtf \
    gtf.list

for SRA_ACCESSION in SRX3632881 SRX3632879 SRX3632877; do

    stringtie -eB \
        $SRA_ACCESSION"_sorted.bam" \
        -o $(pwd)/$SRA_ACCESSION/$SRA_ACCESSION"_schHaem_stringtie-merged.gtf" \
        -A $(pwd)/$SRA_ACCESSION/$SRA_ACCESSION"_schHaem_stringtie-merged_gene_abund.tab" \
        -G schHaem_stringtie-merged.gtf

done

#manually identify the Smp_127030 paralogue
cat ../../dn_ds/Smp_127030_cds_schHae.bed
#KL250964.1      34303   34549   Smp_127030.1.cds
#<...>
#KL250964.1      67553   67895   Smp_127030.1.cds

#Smp_127030 goes from KL250964.1:34303-67895

grep KL250964 SRX3632881/SRX3632881_schHaem_stringtie-merged_gene_abund.tab 
#MSTRG.12886     -       KL250964.1      -       960     31802   18.656673       14.235166       15.798444
#MSTRG.12887     -       KL250964.1      +       34713   65600   0.000000        0.000000        0.000000
#MSTRG.12888     -       KL250964.1      -       38255   40445   0.000000        0.000000        0.000000
#MSTRG.12889     -       KL250964.1      .       60485   60856   0.000000        0.000000        0.000000
#MSTRG.12890     -       KL250964.1      +       66510   67679   0.000000        0.000000        0.000000
#MSTRG.12891     -       KL250964.1      .       78480   78698   64.191780       42.928440       47.642757
#MSTRG.12892     -       KL250964.1      +       101331  114330  18.730738       23.561680       26.149178
#MSTRG.12893     -       KL250964.1      -       113808  114025  0.000000        0.000000        0.000000
#MSTRG.12894     -       KL250964.1      -       126959  141558  3.747135        2.505908        2.781101
#MSTRG.12895     -       KL250964.1      +       142954  143615  5.988595        4.004890        4.444700
#MSTRG.12896     -       KL250964.1      -       148885  149753  35.174911       23.523325       26.106611
#MSTRG.12897     -       KL250964.1      -       149870  152870  19.816387       24.708218       27.421625
#MSTRG.12898     -       KL250964.1      +       162565  178327  85.115547       56.921268       63.172253
#MSTRG.12899     -       KL250964.1      -       176712  200375  15.281089       10.219273       11.341534


#based on these expression values it seems pretty clear that MSTRG.12887 is the
#   right match

grep MSTRG.12887 SRX36328*/*gene_abund.tab
#SRX3632877/SRX3632877_schHaem_stringtie-merged_gene_abund.tab:MSTRG.12887       -       KL250964.1      +       34713   65600   7.318426        3.7638944.831543
#SRX3632879/SRX3632879_schHaem_stringtie-merged_gene_abund.tab:MSTRG.12887       -       KL250964.1      +       34713   65600   3.015467        2.3358052.101585
#SRX3632881/SRX3632881_schHaem_stringtie-merged_gene_abund.tab:MSTRG.12887       -       KL250964.1      +       34713   65600   0.000000        0.0000000.000000


#cleaned up for R
echo "Gene,ID,Gene Name,Reference Strand,Start,End,Coverage,FPKM,TPM"    >schHae_MSTRG12887.csv
echo "MSTRG.12887,-,KL250964.1,+,34713,65600,7.318426,3.7638944.831543" >>schHae_MSTRG12887.csv
echo "MSTRG.12887,-,KL250964.1,+,34713,65600,3.015467,2.3358052.101585" >>schHae_MSTRG12887.csv
echo "MSTRG.12887,-,KL250964.1,+,34713,65600,0.000000,0.0000000.000000" >>schHae_MSTRG12887.csv


#--------------------------- U N U S E D ---------------------------------------
##if we decide to look at expression from each invadolysin paralogue
##get invadolysins (from Berriman 2009 (supp table 8)
#echo "Smp_173070" >berriman_et_al_2009_invadolysin.list
#echo "Smp_090110">>berriman_et_al_2009_invadolysin.list
#echo "Smp_127030">>berriman_et_al_2009_invadolysin.list
#echo "Smp_135530">>berriman_et_al_2009_invadolysin.list
#echo "Smp_090100">>berriman_et_al_2009_invadolysin.list
#echo "Smp_153930">>berriman_et_al_2009_invadolysin.list
#echo "Smp_190480">>berriman_et_al_2009_invadolysin.list
#echo "Smp_167110">>berriman_et_al_2009_invadolysin.list
#echo "Smp_167070">>berriman_et_al_2009_invadolysin.list
#echo "Smp_171340">>berriman_et_al_2009_invadolysin.list
#echo "Smp_167100">>berriman_et_al_2009_invadolysin.list
#echo "Smp_167090">>berriman_et_al_2009_invadolysin.list
#echo "Smp_167120">>berriman_et_al_2009_invadolysin.list
#echo "Smp_171330">>berriman_et_al_2009_invadolysin.list




