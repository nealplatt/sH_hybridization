#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# set_env.sh - sets environmental variables for all remaining scripts.

# TODO(nplatt):

# SET UP MAJOR DIRECTORIES -----------------------------------------------------
HOME_DIR="/master/nplatt/sH_hybridization"

SCRIPTS_DIR="$HOME_DIR/scripts"
DATA_DIR="$HOME_DIR/data"
RESULTS_DIR="$HOME_DIR/results"

RAW_READS_DIR="$DATA_DIR/raw_exome_seq_reads"
GENOME_DIR="$DATA_DIR/genome"

FILTER_READS_DIR="$RESULTS_DIR/filter_reads"
MAP_DIR="$RESULTS_DIR/map_reads"
INTERVALS_DIR="$RESULTS_DIR/intervals"
BQSR_DIR="$RESULTS_DIR/base_recalibration"

SUB_SCRIPTS_DIR="$RESULTS_DIR/sub_scripts"
LOGS_DIR="$RESULTS_DIR/logs"

# SET UP FREQ USED VARS --------------------------------------------------------
REFERENCE="$GENOME_DIR/schHae_v1.fa"
QSUB="qsub -V -cwd -j y -S /bin/bash -q all.q"

# SINGULARITY  -----------------------------------------------------------------
SINGULARITY_IMG="$HOME_DIR/snpCalling_v0.0.6.img"
SINGULARITY="singularity exec $SINGULARITY_IMG"

# CREATE DOWNSTREAM LISTS ------------------------------------------------------
SAMPLE_LIST="$RESULTS_DIR/sample.list"
ls $RAW_READS_DIR >$SAMPLE_LIST

# DELETE LOG AND SCRIPT --------------------------------------------------------
DELETE () {
    #delete log file    
    if [ -f "$1" ]; then
        rm $1
    fi
    
    #delete script file
    if [ -f "$2" ]; then
        rm $2
    fi
}

SUBMIT () {
    echo $1 >$2
    cat $2 | $3
}


WAIT_FOR_CLEAR_QUEUE () {
    NUM_JOBS_IN_QUEUE=$(qstat | grep "snp." | wc -l)
   
    while [ $NUM_JOBS_IN_QUEUE -gt 0 ]; do
        sleep 60s
        echo -n "."
        NUM_JOBS_IN_QUEUE=$(qstat | grep $1 | wc -l)
    done
}

LIMIT_RUNNING_JOBS_TO () {
    # $1 = num max jobs
    NUM_JOBS_IN_QUEUE=$(qstat | grep "snp." | wc -l)
   
    while [ $NUM_JOBS_IN_QUEUE -gt $1 ]; do
        sleep 60s
        echo -n "."
        NUM_JOBS_IN_QUEUE=$(qstat | grep "snp." | wc -l)
    done
}
