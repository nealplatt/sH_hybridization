#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# get_sra_data.sh - download data from SRA, in addition copies over data from 
#    fred's work to use as mansoni outgroup samples

# TODO(neal):

source /master/nplatt/sH_hybridization/scripts/set-env.sh

mkdir $SRA_DIR/
cd $SRA_DIR

#bovis
#ERR119622
#ERR103048	
#ERR119622

#currasoni
#ERR310937 - failing for some reason
#ERR119623

#mansoni
#Sm.BR_0447.1
#Sm.BR_2039.1
#Sm.BR_1278.

#haem
#ERR084970	
#ERR037800	
#SRR433865

for SRA in ERR119622 ERR103048 ERR119622 ERR310937 ERR119623 ERR084970 ERR037800 SRR433865; do
    $SINGULARITY fastq-dump --split-files --gzip $SRA &
done

#get fred's exome data
cp /data/infectious/schistosome/04\ -\ Population\ genetics\ \(Sm\ -\ Sro\)/2017-01-20\ Make\ Analysis\ Global\ Again/v7/1-Alignment/data/Sm.BR_0447.1/*.fastq.gz . 
cp /data/infectious/schistosome/04\ -\ Population\ genetics\ \(Sm\ -\ Sro\)/2017-01-20\ Make\ Analysis\ Global\ Again/v7/1-Alignment/data/Sm.BR_2039.1/*.fastq.gz . 
cp /data/infectious/schistosome/04\ -\ Population\ genetics\ \(Sm\ -\ Sro\)/2017-01-20\ Make\ Analysis\ Global\ Again/v7/1-Alignment/data/Sm.BR_1278.1/*.fastq.gz . 

wait 

rename _1.fastq.gz _R1.fastq.gz *_1.fastq.gz
rename _2.fastq.gz _R2.fastq.gz *_2.fastq.gzls -lh
