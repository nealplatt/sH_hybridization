setwd("C:/Users/Platts-home/Dropbox/work/projects/sH_hybridization/results/bigFig/")

################## D
d<-read.csv("d_2018-06-12.csv", header=T)

d_cols=c(
  rep("black", 3860),
  rep("red", 1851),
  rep("black", 2609),
  rep("red", 2669),
  rep("black", 1400),
  rep("red", 1367),
  rep("black", 667)
)

plot(d$cul_pos, d$D100, ylim=c(-1,1), col=d_cols, pch=19, ylab="D")


#################### FST

#SM_V7_1	3049
#SM_V7_2	1571
#SM_V7_3	1582
#SM_V7_4	1580
#SM_V7_5	683
#SM_V7_6	867
#SM_V7_7	624

fst<-read.csv("fst_2018-06-12.csv", header=T)
fst_cols=c(
  rep("black", 3049),
  rep("red", 1571),
  rep("black", 1582),
  rep("red", 1580),
  rep("black", 683),
  rep("red", 867),
  rep("black", 624)
)

plot(fst$CUL_POS, fst$WEIGHTED_FST, ylim=c(0,1), col=fst_cols, pch=19, ylab="FST")
plot(fst$CUL_POS, fst$MEAN_FST, ylim=c(0,1), col=fst_cols, pch=19, ylab="FST")

#################### SPRIME


s<-read.csv("sprime_2018-06-12.csv", header=T)
s_cols=c(
  rep("black", 77),
  rep("red", 35),
  rep("black", 38),
  rep("red", 29),
  rep("black", 5),
  rep("red", 9),
  rep("black", 11)
)

plot(s$CUL_START, s$SCORE, col=s_cols, ylim=c(0,750000), pch=19, ylab="S")


#################### PCAdmix


ad<-read.csv("pcadmix_MIN_2018-06-12.csv", header=T)
ad_cols=c(
  rep("black", 1547),
  rep("red", 758),
  rep("black", 763),
  rep("red", 667),
  rep("black", 254),
  rep("red", 310),
  rep("black", 279)
)

plot(ad$CUL_POS,ad$P_BOV, col=ad_cols, ylim=c(0,1), pch=19, ylab="Percent S. bovis alleles")




