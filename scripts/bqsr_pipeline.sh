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
ROUND=r1
IN_BAM=$MAP_DIR/$SAMPLE"_processed.bam"
PREV_ROUND_COV_TABLES=" NEED TO FIGURE OUT "


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
