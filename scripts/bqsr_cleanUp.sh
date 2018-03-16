#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr_cleanUp.sh - remove unneccessary files from this round of bqsr

# TODO(nplatt): verify dirs and files after first pass.

cd $BQSR_DIR

# CLEANUP ----------------------------------------------------------------------
trash-put $ROUND"_hc"
trash-put $ROUND"_individual_vcf"
trash-put $ROUND"_db"
trash-put $ROUND"_genotype"
trash-put $ROUND"_cohort_vcf"

#compress and archive scripts and logs
