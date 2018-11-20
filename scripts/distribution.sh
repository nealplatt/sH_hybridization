 grep  "Sh.NE"  ../dating/blocks.tsv \
    | awk '{print $2"\t"$3"\t"$4}' \
    | sed 1d \
    | bedtools sort \
    >bovis_blocks_per_indiv.bed

bedtools merge \
    -i bovis_blocks_per_indiv.bed \
    >bovis_blocks.bed


bedtools merge -d 25000 -i ../tmrca/schHae_v1_probes_schMan_lifted_sorted.bed >merged_probes_25k.bed                                                                                 
bedtools intersect -c -a 1mW_1mS_schMan.bed -b merged_probes_25k.bed >windows_v_probes.bed                                                                                           
awk '{if ($4<5) print $0}' windows_v_probes.bed >1mW_1mS_windows_v_probes_gt5_25k.bed

#manually removed unused contigs (haplotypes ZW)

#count the number of bovis "blocks"
bedtools intersect -c -a 1mW_1mS_windows_v_probes_gt5_25k.bed -b bovis_blocks_per_indiv.bed >introgression_in_windows.bed

#find regions lacking bovis blocks
cat introgression_in_windows.bed | awk '{if ($5==0) print $0}' >no_int_windows.bed

#merge neighboring regions with out bovis
bedtools merge -d 1 -i no_int_windows.bed >merged.bed

#calculate sizes of regions
awk '{print $0"\t"$3-$2}' merged.bed >merged_sizes.bed



#sort
sort -n -k4 merged_sizes.bed

wc -l 1mW_1mS_windows_v_probes_gt5_25k.bed


#make genome file
../../data/genome/schMan_v7.fa.fai
#randomize
awk '{if ($4<5) print $0}' windows_v_probes.bed >1mW_1mS_windows_v_probes_lt5_25k.bed


for REP in $(seq -w 1 1000); do
   
 bedtools shuffle \
        -i bovis_blocks.bed \
        -g schMan_v7.genome \
        -excl 1mW_1mS_windows_v_probes_lt5_25k.bed \
        -seed $RANDOM \
        -noOverlapping \
         >random_blocks.bed

    bedtools intersect \
        -c \
        -a 1mW_1mS_windows_v_probes_gt5_25k.bed \
        -b random_blocks.bed \
        >random_in_windows.bed


    cat random_in_windows.bed | awk '{if ($5==0) print $0}' >random_no_int_windows.bed

    bedtools merge -i random_no_int_windows.bed | awk '{print $0"\t"$3-$2}' >random_merged.bed

    MAX=$(cat random_merged.bed | sort -nr -k4,4 | head -n1 | cut -f4)
    MD5=$(md5sum random_merged.bed | cut -f1 -d" ")
    echo -e $REP"\t"$MAX"\t"$MD5

    cat random_merged.bed >>all.bed

    rm random_blocks.bed
    rm random_in_windows.bed
    rm random_no_int_windows.bed
    rm random_merged.bed

done >max_sizes_random.tab

#to see max
sort -nr -k2 max_sizes_random.tab | head


