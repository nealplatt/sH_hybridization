#clean and process reads to the haematobium genome
source activate snp_calling
source /master/nplatt/schisto_hybridization/scripts/set_env.sh


cd $RESULTS_DIR
mkdir mansoni_trans
cd mansoni_trans

SAMPLES=(   "ERR022875"     "ERR022876"     "ERR022877"     "ERR022878"     
            "ERR022879"     "ERR022880"     "ERR022881"     "ERR022882"
            "ERR022883"     "ERR022874"     "ERR022873"     "ERR022872" )

mkdir logs

#get the data
for SAMPLE in "${SAMPLES[@]}"; do
    fastq-dump --split-files --outdir $SAMPLE"_reads" --gzip $SAMPLE
done

#SAMPLE INFO
#===============================================================================
#Run	    Experiment	InsertSize	MBases	TimeUnit	age	    dev_stage      # 
#-------------------------------------------------------------------------------
#ERR022883	somule6	    490	        3,549	hours	    24	    larva
#ERR022882	somule5	    466	        3,225	hours	    24	    larva
#ERR022881	somule4	    471	        3,655	hours	    24	    larva
#ERR022880	somule3	    467	        3,668	hours	    24	    larva
#ERR022879	somule2	    474	        4,249	hours	    3	    larva
#ERR022878	cerc13	    350	        3,170	hours	    4	    larva
#ERR022877	cerc12	    350	        4,461	hours	    4	    larva
#ERR022876	somule1	    350	        3,439	hours	    3	    larva
#ERR022875	tail1	    350	        1,884	hours	    0	    larva
#ERR022874	somule1	    350	        1,021	hours	    3	    larva
#ERR022873	adult2	    350	        1,525	weeks	    7	    adult
#ERR022872	cerc10a	    175	        2,442	hours	    4	    larva
#===============================================================================

#map reads to the S. mansoni (v7) genome
GENOME=/master/nplatt/schisto_hybridization/data/genome/schMan_v7.fa

for SAMPLE in "${SAMPLES[@]}"; do

     tophat \
         --mate-inner-dist 350 \
         --num-threads 12 \
         --GTF /master/nplatt/schisto_hybridization/data/Sm_v7.0.gff \
         --output-dir $SAMPLE"_tophat" \
         schMan_v7 \
         $SAMPLE"_reads"/$SAMPLE"_1.fastq.gz" \
         $SAMPLE"_reads"/$SAMPLE"_2.fastq.gz" \
         >>logs/$SAMPLE"_tophat.log" 2>&1 &
 
     wait
done

#cufflinks
for SAMPLE in "${SAMPLES[@]}"; do

    cufflinks \
        --no-update-check \
        --num-threads 6 \
        --GTF /master/nplatt/schisto_hybridization/data/Sm_v7.0.gff \
        --multi-read-correct \
        --output-dir $SAMPLE"_cufflinks" \
        --frag-bias-correct $GENOME \
        --library-type fr-unstranded \
        --upper-quartile-norm \
        $SAMPLE"_tophat"/accepted_hits.bam \
         >>$SAMPLE"_cufflinks.log" 2>&1 &

      wait
done

ls ERR*cufflinks/transcripts.gtf >assembly_GTF_list.txt

cuffmerge \
    -o cuffmerge \
    --num-threads 12 \
    --ref-gtf /master/nplatt/schisto_hybridization/data/Sm_v7.0.gff \
    --ref-sequence $GENOME \
    assembly_GTF_list.txt \
    >cuffmerge.log 2>&1 &

tail -f cuffmerge.log

cuffquant \
    --no-update-check \
    --output-dir cuffquant \
    --num-threads 12 \
    --library-type fr-unstranded \
    --frag-bias-correct $GENOME \
    --multi-read-correct \
    ./cuffmerge/merged.gtf \
    ./ERR022872_tophat/accepted_hits.bam \
    ./ERR022873_tophat/accepted_hits.bam \
    ./ERR022874_tophat/accepted_hits.bam \
    ./ERR022875_tophat/accepted_hits.bam \
    ./ERR022876_tophat/accepted_hits.bam \
    ./ERR022877_tophat/accepted_hits.bam \
    ./ERR022878_tophat/accepted_hits.bam \
    ./ERR022879_tophat/accepted_hits.bam \
    ./ERR022880_tophat/accepted_hits.bam \
    ./ERR022881_tophat/accepted_hits.bam \
    ./ERR022882_tophat/accepted_hits.bam \
    ./ERR022883_tophat/accepted_hits.bam \
    >cuffquant.log 2>&1 &

tail -f cuffquant.log


cuffnorm \
    --no-update-check \
    --output-dir cuffnorm_Sm_V7 \
    --num-threads 12 \
    --library-type fr-unstranded \
    --compatible-hits-norm \
    --output-format simple-table \
    --library-norm-method classic-fpkm \
    --labels cerc10a-4h,adult2-7w,somule1-3h,tail1-0h,somule1-3h,cerc12-4h,cerc13-4h,somule2-3h,somule3-24h,sommule4-24h,somule5-24h,somule6-24h \
    /master/nplatt/schisto_hybridization/data/Sm_v7.0.gff  \
    ./ERR022872_tophat/accepted_hits.bam \
    ./ERR022873_tophat/accepted_hits.bam \
    ./ERR022874_tophat/accepted_hits.bam \
    ./ERR022875_tophat/accepted_hits.bam \
    ./ERR022876_tophat/accepted_hits.bam \
    ./ERR022877_tophat/accepted_hits.bam \
    ./ERR022878_tophat/accepted_hits.bam \
    ./ERR022879_tophat/accepted_hits.bam \
    ./ERR022880_tophat/accepted_hits.bam \
    ./ERR022881_tophat/accepted_hits.bam \
    ./ERR022882_tophat/accepted_hits.bam \
    ./ERR022883_tophat/accepted_hits.bam \
    >cuffnorm_Sm_V7.log 2>&1 &

# M8/Invadolysin genes of interest
#----------------------------------
# Smp_173070
# Smp_090110
# Smp_127030
# Smp_135530
# Smp_090100
# Smp_153930
# Smp_190480
# Smp_167110
# Smp_167070
# Smp_171340
# Smp_167100
# Smp_167090
# Smp_167120
# Smp_171330


#hitting memory issues trying with fewer samples
#ERR022873	adult2	    350	        1,525	weeks	    7	    adult
#ERR022874	somule1	    350	        1,021	hours	    3	    larva
#ERR022878	cerc13	    350	        3,170	hours	    4	    larva
#ERR022883	somule6	    490	        3,549	hours	    24	    larva

cuffnorm \
    --no-update-check \
    --output-dir cuffnorm_Sm_V7_reducedN \
    --num-threads 12 \
    --library-type fr-unstranded \
    --compatible-hits-norm \
    --output-format simple-table \
    --library-norm-method classic-fpkm \
    --labels adult2,somule1,cerc13,somule6 \
    /master/nplatt/schisto_hybridization/data/Sm_v7.0.gff  \
    ./ERR022873_tophat/accepted_hits.bam \
    ./ERR022874_tophat/accepted_hits.bam \
    ./ERR022878_tophat/accepted_hits.bam \
    ./ERR022883_tophat/accepted_hits.bam \
    >cuffnorm_Sm_V7_reducedN.log 2>&1 &

echo "Smp_173070" >invadolysin_ids.txt
echo "Smp_090110">>invadolysin_ids.txt
echo "Smp_127030">>invadolysin_ids.txt
echo "Smp_135530">>invadolysin_ids.txt
echo "Smp_090100">>invadolysin_ids.txt
echo "Smp_153930">>invadolysin_ids.txt
echo "Smp_190480">>invadolysin_ids.txt
echo "Smp_167110">>invadolysin_ids.txt
echo "Smp_167070">>invadolysin_ids.txt
echo "Smp_171340">>invadolysin_ids.txt
echo "Smp_167100">>invadolysin_ids.txt
echo "Smp_167090">>invadolysin_ids.txt
echo "Smp_167120">>invadolysin_ids.txt
echo "Smp_171330">>invadolysin_ids.txt

grep -f invadolysin_ids.txt cuffnorm_Sm_V7_reducedN/genes.fpkm_table >invadolysins_gene_fpkm.table



#SAMPLE INFO
#===============================================================================
#Run	    Experiment	InsertSize	MBases	TimeUnit	age	    dev_stage      # 
#-------------------------------------------------------------------------------
#ERR022883	somule6	    490	        3,549	hours	    24	    larva
#ERR022882	somule5	    466	        3,225	hours	    24	    larva
#ERR022881	somule4	    471	        3,655	hours	    24	    larva
#ERR022880	somule3	    467	        3,668	hours	    24	    larva
#ERR022879	somule2	    474	        4,249	hours	    3	    larva
#ERR022878	cerc13	    350	        3,170	hours	    4	    larva
#ERR022877	cerc12	    350	        4,461	hours	    4	    larva
#ERR022876	somule1	    350	        3,439	hours	    3	    larva
#ERR022875	tail1	    350	        1,884	hours	    0	    larva
#ERR022874	somule1	    350	        1,021	hours	    3	    larva
#ERR022873	adult2	    350	        1,525	weeks	    7	    adult
#ERR022872	cerc10a	    175	        2,442	hours	    4	    larva
#===============================================================================
 grep -f invadolysin_ids.txt /master/nplatt/schisto_hybridization/data/Sm_v7.0.gff >invadolysins_Sm_v7.0.gff


cuffnorm \
    --no-update-check \
    --output-dir cuffnorm_Sm_V7_only_invadolysins \
    --num-threads 12 \
    --library-type fr-unstranded \
    --compatible-hits-norm \
    --output-format simple-table \
    --library-norm-method classic-fpkm \
    --labels cerc10a-4h,adult2-7w,somule1-3h,tail1-0h,somule1-3h,cerc12-4h,cerc13-4h,somule2-3h,somule3-24h,sommule4-24h,somule5-24h,somule6-24h \
    ./invadolysins_Sm_v7.0.gff \
    ./ERR022872_tophat/accepted_hits.bam \
    ./ERR022873_tophat/accepted_hits.bam \
    ./ERR022874_tophat/accepted_hits.bam \
    ./ERR022875_tophat/accepted_hits.bam \
    ./ERR022876_tophat/accepted_hits.bam \
    ./ERR022877_tophat/accepted_hits.bam \
    ./ERR022878_tophat/accepted_hits.bam \
    ./ERR022879_tophat/accepted_hits.bam \
    ./ERR022880_tophat/accepted_hits.bam \
    ./ERR022881_tophat/accepted_hits.bam \
    ./ERR022882_tophat/accepted_hits.bam \
    ./ERR022883_tophat/accepted_hits.bam \
    >cuffnorm_Sm_V7_only_invadolysins.log 2>&1 &

fastq-dump --split-files --outdir SRR922067_reads --gzip SRR922067	


 grep -f invadolysin_ids.txt TPM_isoforms_Sm_v7.1.tsv >fred_invadolysins.tsv

all_invadolysin_ids.txt


grep -f berriman_et_al_2009_invadolysin.list ../../data/Sm_v7.0.gff | grep previous | cut -f9 | cut -f1 -d ";" | sed 's/ID=//' | sort -u
#Smp_247870 = Smp_135530, Smp_177600
#Smp_247860 = Smp_135530
#Smp_247850 = Smp_135530
#Smp_303070 = Smp_173070
#Smp_303760 = Smp_173070
#Smp_314000 = Smp_167120
#Smp_314010 = Smp_167120
#Smp_336180 = Smp_167110
#Smp_334410 = Smp_167090, Smp_171320, Smp_204600
#Smp_331850 = Smp_167100

nano updated_smps.list


cat updated_smps.list berriman_et_al_2009_invadolysin.list | sort -u >invadolysin_smps.list

head -n1 TPM_isoforms_Sm_v7.1.tsv >schMan_invadolysin_TPM.tsv   
grep -f invadolysin_smps.list TPM_isoforms_Sm_v7.1.tsv >>schMan_invadolysin_TPM.tsv                                      


transcript_id   miracidia       sporocyst_48h   cercariae       somule_3h       somule_24h      adult_m_21d     adult_f_21d     adult_m_28d     adult_f_28d     adult_m_38d     adult_f_38d
Smp_090100.1    0.00    0.00    0.57    0.73    0.30    0.04    0.03    0.02    0.02    0.00    0.00
Smp_090110.1    0.01    0.08    0.00    0.00    0.00    0.18    0.10    0.28    0.13    0.24    0.27
Smp_127030.1    0.11    0.50    0.20    1.38    3.48    6.25    5.24    5.89    5.89    1.64    6.09
Smp_153930.1    1.31    2.94    0.03    1.02    0.37    1.21    7.95    18.44   35.67   1.45    60.88
Smp_167070.1    27.50   6.44    2.89    2.26    0.44    5.27    4.61    3.45    3.41    2.36    4.50
Smp_167070.2    46.39   0.79    0.18    0.00    0.00    0.28    0.05    0.22    0.15    0.74    0.33
Smp_247850.1    0.32    0.50    3.08    7.50    1.76    12.10   41.54   29.30   39.58   3.24    50.50
Smp_247860.1    1.79    0.11    12.37   1.08    0.46    0.53    0.42    0.33    0.40    0.33    0.37
Smp_247870.1    25.06   15.78   40.09   34.52   19.44   7.92    7.42    8.25    8.62    2.86    10.44
Smp_303070.1    312.93  269.86  1.04    0.24    0.46    0.00    0.00    0.00    0.00    0.00    0.00
Smp_303760.1    256.64  269.86  1.04    0.24    0.46    0.00    0.00    0.00    0.00    0.00    0.00
Smp_303760.2    7.24    27.52   1.94    2.26    0.00    0.70    0.53    0.34    0.59    0.22    0.56
Smp_314000.1    1.51    0.18    0.26    0.17    0.75    0.00    0.03    0.05    0.00    0.05    0.00
Smp_314010.1    0.26    0.00    0.00    0.15    0.00    0.06    0.18    0.15    0.15    0.13    0.00
Smp_314010.2    0.02    0.01    0.11    0.05    0.06    0.01    0.00    0.00    0.00    0.00    0.00
Smp_331850.1    9.29    0.78    0.07    0.05    0.00    0.05    0.03    0.00    0.05    0.02    0.00
Smp_334410.1    0.03    0.03    0.16    0.73    0.94    1.71    1.71    0.92    1.25    1.13    0.90
Smp_336180.1    5.06    0.46    0.24    0.00    0.31    0.56    0.00    0.18    0.19    0.44    0.20
Smp_336180.2    0.47    0.20    0.26    0.00    0.07    0.22    0.32    0.27    0.32    0.02    0.44
Smp_336180.3    0.09    0.00    0.00    0.00    1.04    0.89    0.64    0.96    0.24    0.00    0.00
Smp_336180.4    0.00    0.00    0.29    0.22    0.02    0.59    0.00    0.00    0.58    0.00    0.34
Smp_336180.5    0.00    0.09    1.22    2.17    6.48    1.83    1.29    0.94    1.42    0.67    1.72

#cleaned up
