for BAM in $(ls ../processed_reads/sh_mapped_reads/*_processed.bam); do
    SAMPLE_ID=$(basename $BAM _processed.bam)

    echo $SAMPLE_ID
    bedtools coverage \
        -a ../../data/schHae_v1_probes.bed \
        -b $BAM \
        -mean \
        >$SAMPLE_ID".exon_cov"

done

for i in $(ls *_processed.stats); do
    SAMPLE=$(basename $i _processed.stats)
    READS=$(cat $i | grep mapped | head -n1 | cut -f1 -d" ")
    echo -e $SAMPLE"\t"$READS
done





#fyi 



  | bedtools coverage -a ../../data/schHae_v1_probes.bed -b - -mean >ERR037800.exon_cov &
  | bedtools coverage -a ../../data/schHae_v1_probes.bed -b - -mean >ERR119612.exon_cov &
  | bedtools coverage -a ../../data/schHae_v1_probes.bed -b - -mean >ERR119613.exon_cov &
  | bedtools coverage -a ../../data/schHae_v1_probes.bed -b - -mean >ERR310937.exon_cov &

bedtools intersect -bed  -a ../processed_reads/sh_mapped_reads/ERR037800_processed.bam -b ../../data/schHae_v1_probes.bed >ERR037800.bed &
bedtools intersect -bed  -a ../processed_reads/sh_mapped_reads/ERR119612_processed.bam -b ../../data/schHae_v1_probes.bed >ERR119612.bed &
bedtools intersect -bed  -a ../processed_reads/sh_mapped_reads/ERR119613_processed.bam -b ../../data/schHae_v1_probes.bed >ERR119613.bed &
bedtools intersect -bed  -a ../processed_reads/sh_mapped_reads/ERR310937_processed.bam -b ../../data/schHae_v1_probes.bed >ERR310937.bed &


bedtools coverage -a ../../data/schHae_v1_probes.bed -b ERR037800.bed -mean >ERR037800.exon_cov &
bedtools coverage -a ../../data/schHae_v1_probes.bed -b ERR119612.bed -mean >ERR119612.exon_cov &
bedtools coverage -a ../../data/schHae_v1_probes.bed -b ERR119613.bed -mean >ERR119613.exon_cov &
bedtools coverage -a ../../data/schHae_v1_probes.bed -b ERR310937.bed -mean >ERR310937.exon_cov &


#in R (get mean and num baits with more than 2x coverage)
exon_files<-list.files(".", pattern="*exon_cov")

num_probes<-52237

samples=list(0)
means=list(0)
gt_5s=list(0)

for (i in 1:length(exon_files)){
    cov_stats<-read.table(file=exon_files[i], header=FALSE, sep="\t")

    mean<-round(as.numeric(colMeans(cov_stats[4], na.rm=TRUE)), 2)
    gt_5<-round(sum(cov_stats[4] > 5)/num_probes,2)

    samples<-c(samples, exon_files[i])
    means<-c(means, mean)
    gt_5s<-c(gt_5s, gt_5)
 
    write(c(exon_files[i], mean, gt_5), file="exon_summary", append=TRUE)
}


