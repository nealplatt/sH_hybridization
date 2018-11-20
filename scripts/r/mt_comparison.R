setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/vs_mito")
install.packages("scales")
library("scales")
library("RColorBrewer")
counts=read.csv("counts.csv", header=FALSE)

colors=c(rep("deeppink", 48), rep("blue", 47))

par(mar=c(2,2,2,2))
par(mfrow=c(2,1))
plot(log(counts$V2), log(counts$V3), col=colors, pch=19, xlab="map to bovis", ylab="map to haem", xlim=c(2,12), ylim=c(2,12))
plot(counts$V2, counts$V3, col=alpha(colors, 0.4), pch=19, xlab="map to bovis", ylab="map to haem", xlim=c(0,50000), ylim=c(0,50000))
abline(0,1)

x=