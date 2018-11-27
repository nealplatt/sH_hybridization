#invadolysin_gene_tree

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


