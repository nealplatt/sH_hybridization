library(scales)

setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/snp_density")

counts<-read.table("sman_100k-100k.counts", header=F)
colnames(counts)<-c("chr", "pos_start", "pos_end", "snps", "cul_pos")

chr_cols=c(
  rep("gray81", 889),
  rep("gray39", 482),
  rep("gray81", 501),
  rep("gray39", 473),
  rep("gray81", 253),
  rep("gray39", 250),
  rep("gray81", 193)
  )

mean<-mean(counts$snps)
#2.3/10kb


svg("sman_100k-100k.counts.svg")
  plot(counts$cul_pos, counts$snps, 
       col="black", 
       pch=19, 
       cex=0.75,
       xlab="Chromosome",
       #xlim=c(zoom_start, zoom_end),
       ylab="SNP density (per 100kb)",
       xaxt="n",
       ylim=c(0,250),
       type="l"
       )
  abline(h=mean, col="grey")
  
  
  axis(side=1, 
       at = c(0, 888888981357, 137011725, 187470224, 234750005, 260006124, 284995207, 304283228),
       #labels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "")
       labels=c(rep("", 8))
       )
  abline(h=mean(counts$snps), lty=2, col="black")
dev.off()

