#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr_cleanUp.sh - remove unneccessary files from this round of bqsr

# TODO(nplatt): verify dirs and files after first pass.

cd $BQSR_DIR

# CLEANUP ----------------------------------------------------------------------
rm -r $ROUND"_hc"
rm -r $ROUND"_individual_vcf"
rm -r $ROUND"_db"
rm -r $ROUND"_genotype"
rm -r $ROUND"_cohort_vcf"
mkdir $ROUND"_vcfs"
mkdir $ROUND"_bams"
mkdir $ROUND"_tables"


