setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/lroh")
library(scales)



#read in data
ne<-read.table("ne_roh_regions.bed", header=FALSE, sep = "\t")
tz<-read.table("tz_roh_regions.bed", header=FALSE, sep = "\t")

colnames(ne)<-c("CHROM",
                "AUTO_START",
                "AUTO_END",
                "MIN_START",
                "MAX_END",
                "N_VARIANTS_BETWEEN_MAX_BOUNDARIES",
                "N_MISMATCHES",
                "INDV")

colnames(tz)<-c("CHROM",
                "AUTO_START",
                "AUTO_END",
                "MIN_START",
                "MAX_END",
                "N_VARIANTS_BETWEEN_MAX_BOUNDARIES",
                "N_MISMATCHES",
                "INDV")


proportion<-length(tz$N_VARIANTS_BETWEEN_MAX_BOUNDARIES)/length(ne$N_VARIANTS_BETWEEN_MAX_BOUNDARIES)

svg("roh_snp_counts_2018-07-18.svg", width=4, height=6)
boxplot(tz$N_VARIANTS_BETWEEN_MAX_BOUNDARIES, ne$N_VARIANTS_BETWEEN_MAX_BOUNDARIES, 
        names =c("TZ", "Niger"), 
        col = c("blue", "deeppink"), 
        pch=19, 
        cex=0.5, 
        lty=1,
        ylab="Consequetive homozygous variants",
        width=c(proportion, 1),
        ylim=c(0,4000)
        )
dev.off()

coords<-read.csv("all.coords", header=FALSE, sep = ",")
colnames(coords)<-c("size", "snps")

colors=c(rep("deeppink", 11287),
         rep("blue", 2268))

svg("roh_snp_density_2018-07-18.svg")
plot(coords, 
     pch=19,
     col=colors,
     cex=0.66,
     ylab="Homozygous SNPs",
     xlab="Base pairs"
)

legend(1.2e7, 1000,
       c("Niger", "Tanzania"), 
       col=c("deeppink", "blue"),
       cex = 1,
       lwd = 0, 
       lty = 1,
       pch=19
)

dev.off()
