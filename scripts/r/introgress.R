install.packages("introgress")

library(introgress)
setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/introgress")

admix.gen=read.csv("bovis_zanzibar_fixed_snps_all_introgress-matrix.csv", sep=",", header = FALSE)
loci.data=read.csv("bovis_zanzibar_fixed_snps_all_introgress-loci.csv", sep=",", header = TRUE)


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

triangle.plot(hi.index=h.index, int.het=int.het, pdf=FALSE)
#############################################################################################
admix.gen=read.csv("bovis_haem-sra_fixed_snps_introgress-matrix.csv", sep=",", header=FALSE)
loci.data=read.csv("bovis_haem-sra_fixed_snps_all_introgress-loci.csv", sep=",", header = TRUE)

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

triangle.plot(hi.index=h.index, int.het=int.het, pdf=FALSE)
