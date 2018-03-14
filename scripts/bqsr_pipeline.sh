#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr_pipeline.sh - runs through the bqsr workflow for 1 iteration.
#   from haplotype caller to analyzing covariate tables and generating
#   recalibrated bam files

# TODO(neal): update comments

source /master/nplatt/sH_hybridization/scripts/set-env.sh

cd $BQSR_DIR

#<THESE ARE THE MAJOR VARIABLES THAT NEED TO BE CHANGED FOR EACH ROUND>
export ROUND=r1

#where are the raw/unmodified bams coming from
export UNMODIFIED_BAM_DIR=$MAP_DIR

#ex R1 r1_tables/Sh.TZ_UNG0142.2_PRErecal.r1.table"
#ex R2 r2_tables/Sh.TZ_UNG0142.2_POSTrecal.r2.table"
export PREV_RECAL_TABLE_DIR="$BQSR_DIR/r1_tables/"
export PREV_RECAL_TABLE_EXT="_PRErecal.r1.table"

#make major directories
mkdir $ROUND"_hc"
mkdir $ROUND"_individual_VCF"
mkdir $ROUND"_db"
mkdir $ROUND"_genotype"
mkdir $ROUND"_cohort_VCF"
mkdir $ROUND"_vcfs"
mkdir $ROUND"_bams"
mkdir $ROUND"_tables"

# Workflow
$SCRIPTS_DIR/bqsr_haplotypeCaller.sh
$SCRIPTS_DIR/bqsr_prepForGenotype.sh
#$SCRIPTS_DIR/bqsr_genotype.sh
#$SCRIPTS_DIR/bqsr_buildCohortVcf.sh
#$SCRIPTS_DIR/bqsr_recalibration.sh
