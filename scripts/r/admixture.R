library("RColorBrewer")    

setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/admixture")
Q2=read.table("haem_auto_maf_ld.2.Q")
Q3=read.table("haem_auto_maf_ld.3.Q")
Q4=read.table("haem_auto_maf_ld.4.Q")
Q5=read.table("haem_auto_maf_ld.5.Q")
Q6=read.table("haem_auto_maf_ld.6.Q")
Q7=read.table("haem_auto_maf_ld.7.Q")


names<-scan("haem_auto_maf_ld.list", what = "character" )

bovis<-"green"
haem_sra<-"orange"
haem_tz_ung<-"blue"
haem_niger<-"deeppink"
haem_tz_pem<-"aquamarine"
matt<-"purple"
marg<-"brown"
inter<-"yellow"
cur<-"grey"

svg("haem_auto_maf_admixture.svg")
barplot(t(as.matrix(Q3)), 
        col=c(haem_tz_ung, haem_niger, bovis), 
        ylab="Ancestry", 
        space=0)
dev.off()

