for BAM in $(ls ../processed_reads/sh_mapped_reads/*_processed.bam); do
    SAMPLE_ID=$(basename $BAM _processed.bam)

    echo $SAMPLE_ID
    bedtools coverage \
        -a ../../data/schHae_v1_probes.bed \
        -b $BAM \
        -mean \
        >$SAMPLE_ID".exon_cov"

done

