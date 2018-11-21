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
        --num-threads 6 \
        --multi-read-correct \
        --output-dir $SAMPLE"_cufflinks" \
        --frag-bias-correct schHae_v1.fa \
        --library-type fr-unstranded \
        --upper-quartile-norm \
        $SAMPLE"_tophat"/accepted_hits.bam \
         >>$SAMPLE"_cufflinks.log" 2>&1 &

      wait
done

#----------------------------------------------------
#create a list of all the predicted transcripts
ls ERR*cufflinks/transcripts.gtf >assembly_GTF_list.txt

#then use cuffmerge to merge the predicted transcripts across all samples to
# generate a single/consensus/merged annotation
cuffmerge \
    -o cuffmerge \
    --num-threads 12 \
    --ref-sequence schHae_v1.fa \
    assembly_GTF_list.txt 

#----------------------------------------------------
#quantify expression for each gene across all samples
cuffquant \
    --no-update-check \
    --output-dir cuffquant \
    --num-threads 12 \
    --library-type fr-unstranded \
    --frag-bias-correct schHae_v1.fa \
    --multi-read-correct \
    ./cuffmerge/merged.gtf \
    ./SRX3632881_tophat/accepted_hits.bam \
    ./SRX3632879_tophat/accepted_hits.bam \
    ./SRX3632877_tophat/accepted_hits.bam



#----------------------------------------------------
#normalize expression between samples
cuffnorm \
    --no-update-check \
    --output-dir cuffnorm_schHae_v1 \
    --num-threads 12 \
    --library-type fr-unstranded \
    --compatible-hits-norm \
    --output-format simple-table \
    --library-norm-method classic-fpkm \
    --labels eggs,adult_female,adult_male \
    ########/master/nplatt/schisto_hybridization/data/Sm_v7.0.gff  \
    ./SRX3632881_tophat/accepted_hits.bam \
    ./SRX3632879_tophat/accepted_hits.bam \
    ./SRX3632877_tophat/accepted_hits.bam


#get invadolysins (from Berriman 2009 (supp table 8)
echo "Smp_173070" >berriman_et_al_2009_invadolysin.list
echo "Smp_090110">>berriman_et_al_2009_invadolysin.list
echo "Smp_127030">>berriman_et_al_2009_invadolysin.list
echo "Smp_135530">>berriman_et_al_2009_invadolysin.list
echo "Smp_090100">>berriman_et_al_2009_invadolysin.list
echo "Smp_153930">>berriman_et_al_2009_invadolysin.list
echo "Smp_190480">>berriman_et_al_2009_invadolysin.list
echo "Smp_167110">>berriman_et_al_2009_invadolysin.list
echo "Smp_167070">>berriman_et_al_2009_invadolysin.list
echo "Smp_171340">>berriman_et_al_2009_invadolysin.list
echo "Smp_167100">>berriman_et_al_2009_invadolysin.list
echo "Smp_167090">>berriman_et_al_2009_invadolysin.list
echo "Smp_167120">>berriman_et_al_2009_invadolysin.list
echo "Smp_171330">>berriman_et_al_2009_invadolysin.list

grep -f invadolysin_ids.txt \
    cuffnorm_schHae_v1/genes.fpkm_table \
    >invadolysins_gene_fpkm.table


