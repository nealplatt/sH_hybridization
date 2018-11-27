#gene duplications

#get all platyhelminth m8s plus some outgroup metazoans from pfam:
# outgroups are human, drosophila, celegans, danio

PF01457_selected_sequences.fa

cd /master/nplatt/schisto_hybridization/results/m8_gene_family_expansion

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


