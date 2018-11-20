setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/d")
d<-read.csv("window_D_100snp-50snp_autosomal.csv", header=F)

chr<-fst$V1
cul_pos<-fst$V2
window_pos<-fst$V3
dstat<-fst$V4

chr_cols=c(
  rep("black", 4287),
  rep("grey", 2233),
  rep("black", 2248),
  rep("grey", 2187),
  rep("black", 889),
  rep("grey", 1012),
  rep("black", 680)
)

plot(cul_pos, dstat, pch=19, ylim=c(-1,1), col=chr_cols, xlab = "Chromosome", ylab="D")
