setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/ld")

library(hexbin)

niger_ld<-read.table("niger_SM_V7_1_downsampled_10perc.ld", header=TRUE)
tz_ld<-read.table("tz_SM_V7_1_downsampled_10perc.ld", header=TRUE)


niger_dist<-abs(niger_ld$BP_A - niger_ld$BP_B) 
tz_dist<-abs(tz_ld$BP_A - tz_ld$BP_B) 

png('tz_SM_V7_1_downsampled_10perc.ld.png')
plot(x=tz_dist, y=tz_ld$R2, xlim=c(0,88000000), ylim=c(0,1))
dev.off()

png('niger_SM_V7_1_downsampled_10perc.ld.png')
plot(x=niger_dist, y=niger_ld$R2, xlim=c(0,88000000), ylim=c(0,1))
dev.off()

plot(hexbin(x=tz_chr1_dist, y=tz_chr1_ld$R2, xbins = 10000))

cols <- colorRampPalette(c("darkorchid4","darkblue","green","yellow", "red") )

plot(hexbin(x, y + x*(x+1)/4), main = "Example" ,
     colorcut = seq(0,1,length.out=24),
     colramp = function(n) cols(24) ,
     legend = 0 )




niger_haem_ld<-read.table("niger_KL250487.ld", header=TRUE)
tz_haem_ld<-read.table("tz_KL250487.ld", header=TRUE)

niger_haem_dist<-abs(niger_haem_ld$BP_A - niger_haem_ld$BP_B) 
tz_haem_dist<-abs(tz_haem_ld$BP_A - tz_haem_ld$BP_B) 

png('niger_KL250487.ld.png')
plot(x=tz_haem_dist, y=tz_haem_ld$R2, xlim=c(0,1800000), ylim=c(0,1))
dev.off()

png('tz_KL250487.ld.png')
plot(x=niger_haem_dist, y=niger_haem_ld$R2, xlim=c(0,1800000), ylim=c(0,1))
dev.off()


ld_filtered<-read.table("SM_V7_1_ld-filtered_downsampled_10perc.ld", header=TRUE)

ld_filtered_dist<-abs(ld_filtered$BP_A - ld_filtered$BP_B) 

png('SM_V7_1_ld-filtered_downsampled_10perc.ld.png')
plot(x=ld_filtered_dist, y=ld_filtered$R2, xlim=c(0,88000000), ylim=c(0,1))
dev.off()

niger<-read.table("niger.dist")
hist(niger)



bi<-read.table("tz_KL250487_biallelic.ld", header=TRUE)
bi_dist<-abs(bi$BP_A - bi$BP_B) 
png('tz_KL250487_biallelic.ld.png')
plot(x=bi_dist, y=bi$R2, xlim=c(0,1800000), ylim=c(0,1))
dev.off()


beagle<-read.table("SM_V7_1_beagle_downsampled_5perc.ld", header=TRUE)
beagle_dist<-abs(beagle$BP_A - beagle$BP_B) 
png('SM_V7_1_beagle_downsampled_5perc.ld.png')
plot(x=beagle_dist, y=beagle$R2, xlim=c(0,1800000), ylim=c(0,1))

beagle<-read.table("tz_beagle_SM_V7_1.ld", header=TRUE)
beagle_dist<-abs(beagle$BP_A - beagle$BP_B) 
png('tz_beagle_SM_V7_1_10p_90M_2018-07-10.ld.png')
plot(x=beagle_dist, y=beagle$R2, xlim=c(0,9e7), ylim=c(0,1))
dev.off()

exponential.model <- lm(beagle$R2 ~ beagle_dist)
summary(exponential.model)

distance<-seq(0,2500000, 100)
pred_r2<-exp(predict(exponential.model,list(distance)))

plot(x=beagle_dist, y=beagle$R2, xlim=c(0,3e6), ylim=c(0,1), cex=0.25)
lines(distance, pred_r2,lwd=2, col = "red", xlab = "Time (s)", ylab = "Counts")



setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/ld")

haem_ld<-read.csv("schHae_ld.csv", header=FALSE, sep=",")
man_ld<-read.table("schMan_ld.csv", header=FALSE, sep=",")
chr1_ld<-read.table("1.csv", header=FALSE, sep=",")


plot(hexbin(x=haem_ld$V2, y=haem_ld$V1, xbins = 50))
plot(hexbin(x=man_ld$V2, y=man_ld$V1, xbins = 50))
plot(hexbin(x=chr1_ld$V2, y=chr1_ld$V1, xbins=50))


niger_ld<-read.table("niger.csv", header=FALSE, sep=",")
tz_ld<-read.table("tz.csv", header=FALSE, sep=",")


plot(tz_ld$V2, tz_ld$V1)
plot(niger_ld$V2, niger_ld$V1)

beagle<-read.csv("beagle.csv", header=FALSE, sep=",")
plot(beagle$V2, beagle$V1)


ne<-read.table("beagle_ne.csv", header=FALSE, sep=",")
tz<-read.table("beagle_tz.csv", header=FALSE, sep=",")


plot(tz$V2, tz$V1)
plot(ne$V2, ne$V1)


beagle<-read.table("beagle_tz.ld", header=TRUE)
beagle_dist<-abs(beagle$BP_A - beagle$BP_B) 
png('tz_beagle_SM_V7_1.ld.png')
plot(x=beagle_dist, y=beagle$R2, xlim=c(0,3e6), ylim=c(0,1))
dev.off()



test<-read.table("test.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('test.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,10e7), ylim=c(0,1))
dev.off()

test<-read.table("test_ne.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('test_ne.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,10e7), ylim=c(0,1))
dev.off()


library(dplyr)
library(stringr)
library(ggplot2)



hae<-read.csv("schHae_ld.csv", header=FALSE, sep = ",")
colnames(hae) <- c("rsq","dist")

hae$distc <- cut(hae$dist,breaks=seq(from=min(hae$dist)-1,to=max(hae$dist)+1,by=10000))
hae_avg <- hae %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))

hae_avg <- hae_avg %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

ggplot()+
  geom_point(data=hae_avg,aes(x=start,y=mean),size=0.4,colour="grey20")+
  geom_line(data=hae_avg,aes(x=start,y=mean),size=0.3,alpha=0.5,colour="grey40")+
  labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0,2*10^5,4*10^5,6*10^5,8*10^5,1*10^6,1.2*10^6,1.4*10^6,1.6*10^6),labels=c("0","2e5","4e5","6e5","8e5", "10e5", "12e5", "14e5", "16e5"))+
  theme_bw()


man<-read.csv("schMan_ld.csv", header=FALSE, sep = ",")
colnames(man) <- c("rsq","dist")

man$distc <- cut(man$dist,breaks=seq(from=min(man$dist)-1,to=max(man$dist)+1,by=10000))
man_avg <- man %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))

man_avg <- man_avg %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                              end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                              mid=start+((end-start)/2))

ggplot()+
  geom_point(data=man_avg,aes(x=start,y=mean),size=0.5,colour="grey20")+
  geom_line(data=man_avg,aes(x=start,y=mean),size=0.9,alpha=0.5,colour="grey40")+
  geom_point(data=hae_avg,aes(x=start,y=mean),size=0.5,colour="red1")+
  geom_line(data=hae_avg,aes(x=start,y=mean),size=0.9,alpha=0.5,colour="red3")+
  labs(x="Distance (Bases)",y=expression(LD~(r^{2})))+
  #scale_x_continuous(breaks=c(0,2*10^5,4*10^5,6*10^5,8*10^5,1*10^6,1.2*10^6,1.4*10^6,1.6*10^6),
  #                    labels=c("0","2e5","4e5","6e5","8e5", "10e5", "12e5", "14e5", "16e5"))+
  xlim(0,1.6*10^6)+
  ylim(0,0.6)+
  geom_hline(yintercept=0.2, linetype="dashed", 
             color = "grey20", size=1)+
  theme_bw()

l10<-loess(hae$rsq ~ hae$dist, data=hae, span=0.1)


loess.model<-loess(hae$dist ~ hae$rsq)


plot(hae_avg$mean)


test_dist<-abs(test$BP_A - test$BP_B) 
png('test_ne.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,10e7), ylim=c(0,1))
dev.off()
