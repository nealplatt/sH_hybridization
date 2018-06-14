#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# get_sra_data.sh - download data from SRA, in addition copies over data from 
#    fred's work to use as mansoni outgroup samples

# TODO(neal):

source /master/nplatt/sH_hybridization/scripts/set-env.sh
source activate snp_calling

cd $GENOME_DIR
for SPECIES in BOV CUR GUI INT MAR MAT; do
    wget $LINK &
done

#index all genomes
for GENOME in *.fa; do
    bwa index $GENOME
done

#GET ALL SEQUENCE DATA
cd $SEQ_DIR
for SRA_ACCESSION in ${SRA_SAMPLES[@]}; do
    ${ENVIRONMENTS["SINGULARITY"]} fastq-dump --split-files --gzip $SRA_ACCESSION &
done

wait

#get fred's exome data
cp /data/infectious/schistosome/04\ -\ Population\ genetics\ \(Sm\ -\ Sro\)/2017-01-20\ Make\ Analysis\ Global\ Again/v7/1-Alignment/data/Sm.BR_0447.1/*.fastq.gz . 
cp /data/infectious/schistosome/04\ -\ Population\ genetics\ \(Sm\ -\ Sro\)/2017-01-20\ Make\ Analysis\ Global\ Again/v7/1-Alignment/data/Sm.BR_2039.1/*.fastq.gz . 
cp /data/infectious/schistosome/04\ -\ Population\ genetics\ \(Sm\ -\ Sro\)/2017-01-20\ Make\ Analysis\ Global\ Again/v7/1-Alignment/data/Sm.BR_1278.1/*.fastq.gz . 

wait 

#standardize read names (esp from SRA).
rename _1.fastq.gz _R1.fastq.gz *_1.fastq.gz
rename _2.fastq.gz _R2.fastq.gz *_2.fastq.gz

