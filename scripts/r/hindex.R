install.packages("introgress")

library(introgress)
setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/hindex")

admix.gen=read.csv("matrix.csv", sep=",", header = FALSE)
loci.data=read.csv("loci.csv", sep=",", header = TRUE)

introgress.data<-prepare.data(admix.gen=admix.gen, 
                              loci.data=loci.data, 
                              parental1="H", 
                              parental2="B", 
                              pop.id=TRUE, 
                              ind.id=TRUE, 
                              fixed=TRUE, 
                              sep.rows=FALSE,
                              sep.columns = FALSE)

h.index<-est.h(introgress.data=introgress.data, 
               loci.data=loci.data, 
               ind.touse=NULL,
               fixed=TRUE, 
               p1.allele="H", 
               p2.allele="B")

int.het<-calc.intersp.het(introgress.data=introgress.data)

svg("hindex_2018-07-06.svg", width=7, height=7)
triangle.plot(hi.index=h.index, int.het=int.het, pdf=FALSE)
dev.off()

hist(h.index$h[1:52], ylim=c(0,23), xlim=c(0,0.5), breaks=200, xlab="Hybrid index", col="deeppink")
