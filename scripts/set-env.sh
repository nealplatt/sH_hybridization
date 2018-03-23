#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# set_env.sh - sets environmental variables for all remaining scripts.

# TODO(nplatt):

# SET UP MAJOR DIRECTORIES -----------------------------------------------------
export HOME_DIR="/master/nplatt/sH_hybridization"

export SCRIPTS_DIR="$HOME_DIR/scripts"
export DATA_DIR="$HOME_DIR/data"
export RESULTS_DIR="$HOME_DIR/results"

export RAW_READS_DIR="$DATA_DIR/raw_exome_seq_reads"
export SRA_DIR="$DATA_DIR/sra_reads"
export GENOME_DIR="$DATA_DIR/genome"

export FILTER_READS_DIR="$RESULTS_DIR/filter_reads"
export MAP_DIR="$RESULTS_DIR/map_reads"
export INTERVALS_DIR="$RESULTS_DIR/intervals"
export BQSR_DIR="$RESULTS_DIR/base_recalibration"

export SUB_SCRIPTS_DIR="$RESULTS_DIR/sub_scripts"
export LOGS_DIR="$RESULTS_DIR/logs"

# SET UP FREQ USED VARS --------------------------------------------------------
export REFERENCE="$GENOME_DIR/schHae_v1.fa"
export QSUB="qsub -V -cwd -j y -S /bin/bash -q all.q"

# SINGULARITY  -----------------------------------------------------------------
export SINGULARITY_IMG="$HOME_DIR/snpCalling_v0.0.7.img"
export SINGULARITY="singularity exec $SINGULARITY_IMG"

# CREATE DOWNSTREAM LISTS ------------------------------------------------------
export SAMPLE_LIST="$RESULTS_DIR/sample.list"
ls $RAW_READS_DIR >$SAMPLE_LIST

# DELETE LOG AND SCRIPT --------------------------------------------------------
DELETE () {
    #delete log file    
    if [ -f "$1" ]; then
        trash-put $1
    fi
    
    #delete script file
    if [ -f "$2" ]; then
        trash-put $2
    fi
}

SUBMIT () {
    echo $1 >$2
    cat $2 | $3
}


WAIT_FOR_CLEAR_QUEUE () {
    NUM_JOBS_IN_QUEUE=$(qstat | grep "snp." | awk '$5!="dr"' | wc -l)
   
    while [ $NUM_JOBS_IN_QUEUE -gt 0 ]; do
        sleep 20s
        echo -n "."
        NUM_JOBS_IN_QUEUE=$(qstat | grep "snp." | awk '$5!="dr"' | wc -l)
    done
}

LIMIT_RUNNING_JOBS_TO () {
    # $1 = num max jobs
    NUM_JOBS_IN_QUEUE=$(qstat | grep "snp." | wc -l)
   
    while [ $NUM_JOBS_IN_QUEUE -gt $1 ]; do
        sleep 20s
        echo -n "."
        NUM_JOBS_IN_QUEUE=$(qstat | grep "snp." | wc -l)
    done
}

export -f DELETE
export -f SUBMIT
export -f WAIT_FOR_CLEAR_QUEUE
export -f LIMIT_RUNNING_JOBS_TO
