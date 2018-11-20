library("RColorBrewer")

setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/sfs")

sfs<- read.table("joint_fsfs", header = FALSE, sep = ",")

usfs<-read.csv("usfs_xy.csv", header=FALSE, sep=",", f)

niger_ac<-usfs$V1
tz_ac<-usfs$V2
counts=log(usfs$V3)

rbPal<-colorRampPalette(c('red', 'yellow', 'blue', 'purple', 'pink'))
color<-rbPal(5000)[as.numeric(cut(counts,breaks=5000))]


plot(niger_ac, tz_ac,col=color, pch=15)
log10(usfs$V3)

hist(niger_ac)
hist(tz_ac)
