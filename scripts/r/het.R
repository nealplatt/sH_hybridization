setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results")

het<-read.csv("heterozygosity.csv", header=TRUE, row.names=1)

cols<-c(rep("deeppink", 46), rep("blue", 47))


svg("het.svg", width=3, height=6)
boxplot(het$HET[1:46],het$HET[46:92], 
        col=c("deeppink", "blue"), 
        pch=19, 
        notch = TRUE, 
        names = c("NE", "TZ"),
        ylim=c(0,0.5),
        main="Heterozygosity")
dev.off()

svg("f.svg", width=3, height=6)
boxplot(het$F[1:46],het$F[46:92], 
        col=c("deeppink", "blue"), 
        pch=19, 
        notch = TRUE, 
        names = c("NE", "TZ"),
        ylim=c(-1, 1),
        main="F")
dev.off()



mean(het$F[1:46])
mean(het$F[47:92])

mean(het$HET[1:46])
mean(het$HET[47:92])

