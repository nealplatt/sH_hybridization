library(qqman)

setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results")

#read in the list of introgressed snps from sprime
introgressed_snps<-read.table(file="introgressed_snps.list", header=FALSE)

#read in the fst table
fst<-read.table(file="mainChr_namedSNPS.weir.fst", header=TRUE)
names(fst)<-c("NAMES", "CHROM", "POS", "FST")
fstsubset<-fst[complete.cases(fst),]
SNP<-fstsubset$NAMES
mydf<-data.frame(SNP,fstsubset)

test<-unlist(introgressed_snps)

manhattan(mydf, 
          chr="CHROM", 
          bp="POS", 
          p="FST", 
          snp="SNP", 
          logp=FALSE, 
          ylab="WC Fst", 
          chrlabs =c(1, 2, 3, 4, 5, 6, 7, "ZW"),
          highlight=test
          )

introgressed_snps

#test for normality of fst
library(ggpubr)
ggqqplot(fst$FST)
hist(fst$FST)
