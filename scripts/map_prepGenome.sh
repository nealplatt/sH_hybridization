#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# prep_genome - generates gatk, faidx, and bwa indecies for the genome

# TODO(nplatt): auto-Download the genome from NCBI rather than having it stored.

#load variables
source /master/nplatt/sH_hybridization/scripts/set-env.sh

mkdir $GENOME_DIR $GENOME_DIR/scripts $GENOME_DIR/logs
cd $GENOME_DIR

#download genome from NCBI


#softlink to an more appropriate name
ln -s GCA_000699445.1_SchHae_1.0_genomic.fna $REFERENCE

THREADS=1

#BWA index
JOB_NAME=bwa_index
BWA="$SINGULARITY bwa index $REFERENCE"
BWA_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o logs/$JOB_NAME.log"
echo "$BWA" >scripts/$JOB_NAME.sh
cat scripts/$JOB_NAME.sh | $BWA_QSUB

#GATK
JOB_NAME=gatk_dict
GATK="$SINGULARITY gatk CreateSequenceDictionary -R $REFERENCE -O schHae_v1.dict"
GATK_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o logs/$JOB_NAME.log"
echo "$GATK" >scripts/$JOB_NAME.sh
cat scripts/$JOB_NAME.sh | $GATK_QSUB

#faidx
JOB_NAME=sam_faidx
FAIDX="$SINGULARITY samtools faidx $REFERENCE"
FAIDX_QSUB="$QSUB -pe mpi $THREADS -N $JOB_NAME -o logs/$JOB_NAME.log"
echo "$FAIDX" >scripts/$JOB_NAME.sh
cat scripts/$JOB_NAME.sh | $FAIDX_QSUB


