#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# invadolysin_gene_tree.sh - use the dn_ds invadolsyin fasta file to build a
#   gene tree (all sequences and only uniq haplotypes).

# Uses a conda to manage the enviroment

#Set up the environment
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir invadolysin_gene_tree
cd invadolysin_gene_tree

#get the invadolysin phased haplotypes
cp ../dn_ds/Smp_127030_haplotypes_cds.fasta .
#manually remove empty sequences with nano

#get the schMan invadolysin
cp ../dn_ds/schMan_Smp_127030_combined_cds.fasta .

#combine
cat Smp_127030_haplotypes_cds.fasta schMan_Smp_127030_combined_cds.fasta \
    >Smp_127030.fasta

raxmlHPC \
    -f a \
    -m GTRCAT \
    -p 12345 \
    -x 12345 \
    -s Smp_127030.fasta \
    -# 1000 \
    -n Smp_127030_haplotypes \
    -o schMan

#this tree is a bit unweildy...reduce to uniq haplotypes


#get uniq haplotypes
# get number of times each haplotype is present
cat Smp_127030.fasta \
    | fasta_formatter  -t \
    | cut -f2 \
    | sort \
    | uniq -c \
    | awk '{print $2"\t"$1}' \
    >seq.counts
# get names for each haplotype (from the sequence identifiers)
cat Smp_127030.fasta \
    | fasta_formatter  -t \
    | sort -u -k2,2 \
    | awk '{print $2"\t"$1}' \
    >Smp_127030_uniq_haps.tab

#join those two so that each uniq haplotype has a sample id and the number of
#   times that haplotype is in the dataset
 join seq.counts Smp_127030_uniq_haps.tab  \
    | awk '{print ">"$2"-"$3"\n"$1}' \
    >Smp_127030_uniq_haps.fasta

#build a tree
raxmlHPC-PTHREADS \
    -T 12 \
    -f a \
    -m GTRCAT \
    -p 12345 \
    -x 12345 \
    -s Smp_127030_uniq_haps.fasta \
    -# 1000 \
    -n Smp_127030_uniq_haps \
    -o 1-schMan

#the best tree was cleaned up in inkscape to add colored lines.

