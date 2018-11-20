setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/ld")
library(dplyr)
library(stringr)
library(ggplot2)


#read in data
hae<-read.csv("schHae_ld.csv", header=FALSE, sep = ",")
man<-read.csv("schMan_ld.csv", header=FALSE, sep = ",")

#add column names
colnames(man) <- c("rsq","dist")
colnames(hae) <- c("rsq","dist")

#calculate binned distances and average
hae$distc <- cut(hae$dist,breaks=seq(from=min(hae$dist)-1,to=max(hae$dist)+1,by=1000))
hae_avg <- hae %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
hae_avg <- hae_avg %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                              end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                              mid=start+((end-start)/2))

man$distc <- cut(man$dist,breaks=seq(from=min(man$dist)-1,to=max(man$dist)+1,by=1000))
man_avg <- man %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
man_avg <- man_avg %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                              end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                              mid=start+((end-start)/2))


svg("ld_2017-07-17.svg")
plot(
  y = c(man_avg$median, hae_avg$median), 
  x = c(man_avg$start, hae_avg$start), 
  pch=19, 
  col=c(
    rep("ivory", 86517), 
    rep("gray39", 1351)),
  cex=0.75,
  xlim=c(0,1.6e6),
  ylim=c(0,0.75),
  xlab="Distance (basepairs)",
  ylab="R2")

legend(1.2e6, 0.5,
       c("SM_V7", 
         "SchHae_v1"), 
       col=c("ivory", "gray39"),
       cex = 1,
       lwd = 0, 
       lty = 1,
       pch=19,
       box.lty=2,
       title = "Assembly")
dev.off()


#loess smoothing
man_l66<-loess(man_avg$median ~ man_avg$start, data=man_avg, span=0.66)
hae_l66<-loess(hae_avg$median ~ hae_avg$start, data=hae_avg, span=0.66)


man_l66_smooth<-predict(man_l66, se = TRUE)
hae_l66_smooth<-predict(hae_l66, se = TRUE)

lines(hae_l66_smooth$fit[1:1350], x=hae_avg$start[1:1350], col="blue", lwd=2)
lines(man_l66_smooth$fit, x=man_avg$start, col="red", lwd=2)


lines(hae_avg$start[1:1350], hae_l66_smooth$fit[1:1350] - qt(0.975,hae_l66_smooth$df)*hae_l66_smooth$se, lty=2, col="blue", lwd=2)
lines(hae_avg$start[1:1350], hae_l66_smooth$fit[1:1350] + qt(0.975,hae_l66_smooth$df)*hae_l66_smooth$se, lty=2, col="blue", lwd=2)

lines(man_avg$start[1:153], man_l66_smooth$fit[1:153] - qt(0.975,man_l66_smooth$df)*man_l66_smooth$se, lty=2, col="red", lwd=2)
lines(man_avg$start[1:153], man_l66_smooth$fit[1:153] + qt(0.975,man_l66_smooth$df)*man_l66_smooth$se, lty=2, col="red", lwd=2)

