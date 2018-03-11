#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# check_map_process_bam.sh - checks that all the expected output files were
#   created in the mapping/processing of the bam files

# TODO(nplatt): auto-resubmit rather than generating a list.
# TODO(nplatt): if sample has all of its files...delete intermediate files.

#load variables
source master/nplatt/sH_hybridizationscripts/set_env.sh

cd $MAP_DIR

for SAMPLE in $(cat $SAMPLE_LIST); do 
    #Check to see if these files exist
    for FILE in $MAP_DIR/$SAMPLE"_R1.sai" \
                $MAP_DIR/$SAMPLE"_R2.sai" \
                $MAP_DIR/$SAMPLE"_RX.sai" \
                $MAP_DIR/$SAMPLE"_samSE.bam" \
                $MAP_DIR/$SAMPLE"_samPE.bam" \
                $MAP_DIR/$SAMPLE"_samSE_sorted.bam" \
                $MAP_DIR/$SAMPLE"_samPE_sorted.bam" \
                $MAP_DIR/$SAMPLE"_merged.bam" \
                $MAP_DIR/$SAMPLE"_merged_RGs.bam" \
                $MAP_DIR/$SAMPLE"_merged_RGs_sorted.bam"; do
        if [ ! -f $FILE ]; then
            echo "FAILED $FILE"
        fi
    done
done



