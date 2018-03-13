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

# SET UP FREQ USED VARS --------------------------------------------------------
REFERENCE="$GENOME_DIR/schHae_v1.fa"
QSUB="qsub -V -cwd -j y -S /bin/bash -q all.q"

# SINGULARITY  -----------------------------------------------------------------
SINGULARITY_IMG="$HOME_DIR/snpCalling_v0.0.6.img"
SINGULARITY="singularity exec $SINGULARITY_IMG"

# CREATE DOWNSTREAM LISTS ------------------------------------------------------
SAMPLE_LIST="$RESULTS_DIR/sample.list"
ls $RAW_READS_DIR >$SAMPLE_LIST


# WORKFLOW/PIPELINE ------------------------------------------------------------
# 01] set-env.sh
# 02] map_filterReads.sh
# 03] map_prepGenome.sh
# 04] map_mapAndProcessBam.sh
# 05] map_checkBams.sh
# 06] map_buildListOfIntervals.sh
# 07] bqsr-r1_haplotypeCaller.sh
# 08] bqsr-r1_prepForGenotype.sh
# 09] bqsr-r1_genotype.sh
# 10] bqsr-r1_buildCohortVcf.sh
# 11] bqsr-r1_bqsr.sh
# 12] bqsr-r2_haplotypeCaller.sh
# 13] bqsr-r2_prepForGenotype.sh
# 14] bqsr-r2_genotype.sh
# 15] bqsr-r2_buildCohortVcf.sh
# 16] bqsr-r2_bqsr.sh
# 17] bqsr-r2_buildCohortVcf.sh


#--results
#----filter_reads
#----map_reads
#----intervals
#----base_recalibration
#------logs
#------scripts
#------r0_bqsr_tables
#------r1_bqsr_tables
#------rN_bqsr_tables
#------<modified bams>
