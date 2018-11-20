setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/vs_mito/")

data<-read.table("mito_cov.tab", header=TRUE)

plot(x=data$mito_cov, 
     y=data$filt_reads, 
     pch=19, 
     col=c(
       rep("deeppink", 48), 
       rep("blue", 47)
       ),
     cex=1.2,
     xlab = "Mito. coverage",
     ylab = "Filtered reads (x10^6)"
     )

barplot(data$mito_cov,
        col=c(
          rep("deeppink", 48), 
          rep("blue", 47)),
        ylim=c(0, 225),
        xlab="Samples",
        ylab="Mito. coverage"
        )
