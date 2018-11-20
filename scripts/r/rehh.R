setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/rehh")

install.packages("rehh")

library("rehh")

par(mfrow=c(1,1))


################ 1
ne_hap<-data2haplohh(hap_file=c("niger_1.inp"), 
                     map_file="niger.map",
                     recode.allele=TRUE,
                     chr.name = 1
                    )

tz_hap<-data2haplohh(hap_file="tz_1.inp", 
                     map_file="tz.map",
                     recode.allele=TRUE,
                     chr.name=1
                     )

tz_hh<-scan_hh(tz_hap, threads=2)
ne_hh<-scan_hh(ne_hap, threads=2)

xpehh_1<-ies2xpehh(hh_pop1 = tz_hh,
                   hh_pop2 = ne_hh,
                   popname1="TZ",
                   popname2="NE",
                   method="bilateral")

################ 2
ne_hap<-data2haplohh(hap_file=c("niger_2.inp"), 
                     map_file="niger.map",
                     recode.allele=TRUE,
                     chr.name = 2
)

tz_hap<-data2haplohh(hap_file="tz_2.inp", 
                     map_file="tz.map",
                     recode.allele=TRUE,
                     chr.name=2
)

tz_hh<-scan_hh(tz_hap, threads=2)
ne_hh<-scan_hh(ne_hap, threads=2)

xpehh_2<-ies2xpehh(hh_pop1 = tz_hh,
                   hh_pop2 = ne_hh,
                   popname1="TZ",
                   popname2="NE",
                   method="bilateral")

################ 3
ne_hap<-data2haplohh(hap_file=c("niger_3.inp"), 
                     map_file="niger.map",
                     recode.allele=TRUE,
                     chr.name = 3
)

tz_hap<-data2haplohh(hap_file="tz_3.inp", 
                     map_file="tz.map",
                     recode.allele=TRUE,
                     chr.name=3
)

tz_hh<-scan_hh(tz_hap, threads=2)
ne_hh<-scan_hh(ne_hap, threads=2)

xpehh_3<-ies2xpehh(hh_pop1 = tz_hh,
                   hh_pop2 = ne_hh,
                   popname1="TZ",
                   popname2="NE",
                   method="bilateral")

################ 4
ne_hap<-data2haplohh(hap_file=c("niger_4.inp"), 
                     map_file="niger.map",
                     recode.allele=TRUE,
                     chr.name = 4
)

tz_hap<-data2haplohh(hap_file="tz_4.inp", 
                     map_file="tz.map",
                     recode.allele=TRUE,
                     chr.name=4
)

tz_hh<-scan_hh(tz_hap, threads=2)
ne_hh<-scan_hh(ne_hap, threads=2)

xpehh_4<-ies2xpehh(hh_pop1 = tz_hh,
                   hh_pop2 = ne_hh,
                   popname1="TZ",
                   popname2="NE",
                   method="bilateral")
################ 5
ne_hap<-data2haplohh(hap_file=c("niger_5.inp"), 
                     map_file="niger.map",
                     recode.allele=TRUE,
                     chr.name = 5
)

tz_hap<-data2haplohh(hap_file="tz_5.inp", 
                     map_file="tz.map",
                     recode.allele=TRUE,
                     chr.name=5
)

tz_hh<-scan_hh(tz_hap, threads=2)
ne_hh<-scan_hh(ne_hap, threads=2)

xpehh_5<-ies2xpehh(hh_pop1 = tz_hh,
                   hh_pop2 = ne_hh,
                   popname1="TZ",
                   popname2="NE",
                   method="bilateral")

################ 6
ne_hap<-data2haplohh(hap_file=c("niger_6.inp"), 
                     map_file="niger.map",
                     recode.allele=TRUE,
                     chr.name = 6
)

tz_hap<-data2haplohh(hap_file="tz_6.inp", 
                     map_file="tz.map",
                     recode.allele=TRUE,
                     chr.name=6
)

tz_hh<-scan_hh(tz_hap, threads=2)
ne_hh<-scan_hh(ne_hap, threads=2)

xpehh_6<-ies2xpehh(hh_pop1 = tz_hh,
                   hh_pop2 = ne_hh,
                   popname1="TZ",
                   popname2="NE",
                   method="bilateral")

################ 7
ne_hap<-data2haplohh(hap_file=c("niger_7.inp"), 
                     map_file="niger.map",
                     recode.allele=TRUE,
                     chr.name = 7
)

tz_hap<-data2haplohh(hap_file="tz_7.inp", 
                     map_file="tz.map",
                     recode.allele=TRUE,
                     chr.name=7
)

tz_hh<-scan_hh(tz_hap, threads=2)
ne_hh<-scan_hh(ne_hap, threads=2)

xpehh_7<-ies2xpehh(hh_pop1 = tz_hh,
                   hh_pop2 = ne_hh,
                   popname1="TZ",
                   popname2="NE",
                   method="bilateral")

################
xpehh<-rbind(xpehh_1, xpehh_2, xpehh_3, xpehh_4, xpehh_5, xpehh_6, xpehh_7)


xpehhplot(xpehh)

write.csv(xpehh, file="xpehh_maf00_2018-09-06.csv")


tz_ihs<-ihh2ihs(tz_hh, minmaf=0.05, freqbin=0.01)

ne_ihs<-ihh2ihs(ne_hh, minmaf=0.5, freqbin=0.6)


calc_ehhs(ne_hap, mrk=2425, plotehhs = TRUE, main = "ne")
calc_ehhs(tz_hap, mrk=2425, plotehhs = TRUE, main = "tz")
 
calc_ehhs(ne_hap, mrk=2311, plotehhs = TRUE)
calc_ehhs(tz_hap, mrk=2310, plotehhs = TRUE, )



plot(y=xpehh$`XPEHH (TZ vs. NE)`[18151:19159], x=xpehh$POSITION[18151:19159], pch=19)
