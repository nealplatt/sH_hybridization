
#now get the mansoni sequence and add to use as an outgroup
samtools faidx ../../data/genome/schMan_v7.fa SM_V7_4 >SM_V7_4.fasta

${ENVIRONMENTS["TITAN SINGULARITY"]} \
    /usr/software/progressiveCactus/submodules/hal/bin/halLiftover  \
        $RESULTS_DIR/wga/schMan7_vs_schMan1.hal \
        schHae_v1 \
        Smp_127030_UNprobed-cds_schHae.bed \
        schMan_v7 \
        Smp_127030_UNprobed-cds_schMan.bed

bedtools maskfasta -fi SM_V7_4.fasta -bed Smp_127030_UNprobed-cds_schMan.bed -fo SM_V7_4_masked.fasta

${ENVIRONMENTS["TITAN SINGULARITY"]} \
    /usr/software/progressiveCactus/submodules/hal/bin/halLiftover  \
        $RESULTS_DIR/wga/schMan7_vs_schMan1.hal \
        schHae_v1 \
        Smp_127030_cds_schHae.bed \
        schMan_v7 \
        Smp_127030_cds_schMan.bed

bedtools getfasta \
        -bed Smp_127030_cds_schMan.bed \
        -fi SM_V7_4_masked.fasta \
        -fo schMan_Smp_127030_cds_indiv.fasta

grep -v "^>" schMan_Smp_127030_cds_indiv.fasta \
        | awk 'BEGIN { ORS=""; print ">schMan\n" } { print }' > schMan_Smp_127030_combined_cds.fasta

    cat schMan_Smp_127030_combined_cds.fasta >>Smp_127030_cds.fasta
    echo >>Smp_127030_cds.fasta

