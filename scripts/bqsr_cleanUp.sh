#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bqsr_cleanUp.sh - remove unneccessary files from this round of bqsr

# TODO(nplatt): 

cd $BQSR_DIR

# CLEANUP ----------------------------------------------------------------------

TRASH=$BQSR_DIR/$ROUND"_trash"
# remove intermediate directories
mkdir $TRASH
mkdir $TRASH/$ROUND"_vcfs"

mv $BQSR_DIR/$ROUND"_hc" $TRASH
mv $BQSR_DIR/$ROUND"_individual_vcf" $TRASH
mv $BQSR_DIR/$ROUND"_db" $TRASH
mv $BQSR_DIR/$ROUND"_genotype" $TRASH
mv $BQSR_DIR/$ROUND"_cohort_vcf" $TRASH

#remove intermediate files
mv $BQSR_DIR/$ROUND"_vcfs/cohort_"$ROUND"_rawINDELS.g.vcf"  $TRASH/$ROUND"_vcfs"
mv $BQSR_DIR/$ROUND"_vcfs/cohort_"$ROUND"_rawINDELS.g.vcf.idx" $TRASH/$ROUND"_vcfs"
mv $BQSR_DIR/$ROUND"_vcfs/cohort_"$ROUND"_filteredINDELS.g.vcf" $TRASH/$ROUND"_vcfs"
mv $BQSR_DIR/$ROUND"_vcfs/cohort_"$ROUND"_filteredINDELS.g.vcf.idx" $TRASH/$ROUND"_vcfs"
mv $BQSR_DIR/$ROUND"_vcfs/cohort_"$ROUND"_rawSNPS.g.vcf" $TRASH/$ROUND"_vcfs"
mv $BQSR_DIR/$ROUND"_vcfs/cohort_"$ROUND"_rawSNPS.g.vcf.idx" $TRASH/$ROUND"_vcfs"
mv $BQSR_DIR/$ROUND"_vcfs/cohort_"$ROUND"_filteredSNPS.g.vcf" $TRASH/$ROUND"_vcfs"
mv $BQSR_DIR/$ROUND"_vcfs/cohort_"$ROUND"_filteredSNPS.g.vcf.idx" $TRASH/$ROUND"_vcfs"
mv $BQSR_DIR/$ROUND"_vcfs/merge_variants.list" $TRASH/$ROUND"_vcfs"

#now with all the files organized (safely) delete them
trash-put $TRASH


#compress and archive  logs
cd $LOGS_DIR

ls *.log >log.list

for LOG in $(cat log.list); do
    mv $LOG $ROUND"_logs"
done

tar -cf $ROUND"_logs.tar" $ROUND"_logs"

rm -r log.list $ROUND"_logs"

#compress and archive shell scripts
cd $SUB_SCRIPTS_DIR
mkdir $ROUND"_shellscripts"

ls *.sh >sh.list

for SH in $(cat sh.list); do
    mv $SH $ROUND"_shellscripts"
done

tar -cf $ROUND"_shellscripts.tar" $ROUND"_shellscripts" &
rm -r sh.list $ROUND"_shellscripts"

########################## END CURRENT ROUND OF RECAL ##########################


