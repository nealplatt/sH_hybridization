#mitochondrial assembly

#clean and process reads to the haematobium genome
source activate snp_calling
source /master/nplatt/schisto_hybridization/scripts/set_env.sh

#START FILTERING/MAPPING PROCESS
cd $RESULTS_DIR

FILTER_DIR=$RESULTS_DIR/processed_reads/wgrs
MITO_DIR=$RESULTS_DIR/mitobim

mkdir -p $MITO_DIR

cd $MITO_DIR

#run mitobim (ON TITAN)
for SAMPLE in "${WGRS_SAMPLES[@]}"; do
    #mkdir $SAMPLE

    cd /master/nplatt/schisto_hybridization/results/mitobim/$SAMPLE

    PE_1=/master/nplatt/schisto_hybridization/results/processed_reads/wgrs/$SAMPLE"_wgrs_filtered_R1_paired.fastq.gz"
    PE_2=/master/nplatt/schisto_hybridization/results/processed_reads/wgrs/$SAMPLE"_wgrs_filtered_R2_paired.fastq.gz"

    #zcat $PE_1 $PE_2 >$SAMPLE.fastq &
    #cp ../schHae_mt.fa . &
    rm -r iteration*

    ~/bin/MITObim/MITObim.pl \
        -start 1 \
        -end 100 \
        -sample $SAMPLE \
        -ref schHae \
        -readpool $SAMPLE.fastq  \
        --quick schHae_mt.fa \
        --NFS_warn_only \
        >mitobim_quick.log &
        #| tee mitobim_quick.log

    cd /master/nplatt/schisto_hybridization/results/mitobim
done



#now polish each initial assembly
for SAMPLE in "${WGRS_SAMPLES[@]}"; do
    cd /master/nplatt/schisto_hybridization/results/mitobim/$SAMPLE
    
    INITIAL_ASSEMBLY=$(cat mitobim_quick.log | grep "Final assembly result will be written to file:" | cut -f2 -d":" | sed 's/ //')

    cp $INITIAL_ASSEMBLY .

    INITIAL_ASSEMBLY=$(basename $INITIAL_ASSEMBLY)

    mkdir quick
    mv iteration* quick/

    ~/bin/MITObim/MITObim.pl \
        -start 1 \
        -end 100 \
        -sample $SAMPLE \
        -ref $INITIAL_ASSEMBLY \
        -readpool $SAMPLE.fastq  \
        --quick $INITIAL_ASSEMBLY \
        --NFS_warn_only \
        > mitobim_polish.log &
        #| tee mitobim_polish.log

    cd ..
done




