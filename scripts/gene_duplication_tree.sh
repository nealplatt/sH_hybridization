#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# gene_duplication_tree.sh - build a tree of m8 peptidases to identify
#   gene duplication nodes

# Uses a conda to manage the enviroment

#Set up the environment
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir m8_gene_family_expansion
cd m8_gene_family_expansion

#get all platyhelminth m8s plus some outgroup metazoans from pfam:
# outgroups are human, drosophila, celegans, danio
# DL'd as ->PF01457_selected_sequences.fa

#convert to single line
fasta_formatter -w0 -i PF01457_selected_sequences.fa >SL.fas

#find HExxH
cat SL.fas | grep -B1 HE..H >HExxH.fas

#extract peptidase
python ../../scripts/extract_peptidase.py HExxH.fas HExxH_peptidase.fas

#align
muscle -in HExxH_peptidase.fas -out HExxH_peptidase_muscle.fas

#manually edit headers
#m8_platy_meta_peptidase_HExxH_trimmed.fa

#remove sites with fewer than 75% aligned bases
#m8_platy_meta_75perc.fas

#cleanup for raxml
sed -i 's/:/#/gi' m8_platy_meta_75perc.fas

#build a tree
 raxmlHPC-PTHREADS \
    -T 12 \
    -m PROTGAMMAWAG \
    -p 12345 \
    -s m8_platy_meta_75perc.fas \
    -# 100 \
    -n m8_platy_meta_75perc

#this tree is imported into megas gene duplication wizard and species names maped
# tree is manipulated in raxml and final tree is polished in inkscape


