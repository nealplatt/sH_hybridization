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
SINGULARITY_IMG="$HOME_DIR/snpCalling_v0.0.5.img"
SINGULARITY="singularity exec $SINGULARITY_IMG"

# CREATE DOWNSTREAM LISTS ------------------------------------------------------
SAMPLE_LIST="$RESULTS_DIR/sample.list"
ls $RAW_READS_DIR >$SAMPLE_LIST


# WORKFLOW/PIPELINE ------------------------------------------------------------
# 01] set_env.sh
# 02] filter_reads.sh AND prep_genomes.sh
# 03] map_process_bam.sh
# 04] check_map_process_bam.sh
# 05] build_gatk_intervals.sh
# 06] bsrcl_hc.sh
# 07] bsrcl_to_genotype.sh
# 08] bsrcl_genotype.sh
# 09] build_cohort_vcf.sh
# 10] bsrcl_recal.sh
# 11] bsrcl_hc_r2.sh

