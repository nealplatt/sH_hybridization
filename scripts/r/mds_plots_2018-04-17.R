setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results")

#R code
auto_d <- read.table("auto_mds.mds", h=T)

auto_d$specLoc = factor(c(
                 rep("SRA", 3),
                 rep("NE.Dai", 13),
                 rep("NE.Doki", 1),
                 rep("NE.Kar", 8),
                 rep("NE.Lata", 7),
                 rep("NE.LibTB", 5),
                 rep("NE.NG", 4),
                 rep("NE.Seb", 3),
                 rep("NE.Tiag", 1),
                 rep("NE.YK", 4),
                 rep("NE.Youri", 2),
                 rep("Tz.PEM", 26),
                 rep("Tz.UNG", 21)))

auto_d$pop = factor(c(
                 rep("SRA", 3),
                 rep("NE", 48),
                 rep("TZ_PEM", 26),
                 rep("TZ_UNG", 21)))


plot(auto_d$C1, auto_d$C2, col=c(rep("chartreuse3", 3), rep("deeppink", 48), rep("blue", 26), rep("cyan", 21)), lwd=1, pch=19, cex=1.5, xlab="PC 1", ylab="PC 2", main = "MDS of Sh (gatk/plink)")

legend(0.03, 0.025,
       c("SRA", "Niger", "Tanzania (PEM)", "Tanzania (UNG)"),
       col = c("chartreuse3", "deeppink", "blue", "cyan"),
       cex = 0.8,
       lwd = 7, lty = 1)


mito_d <- read.table("mito_mds.mds", h=T)

mito_d$specLoc = factor(c(
  rep("SRA", 3),
  rep("NE.Dai", 6),
  rep("NE.Kar", 3),
  rep("NE.Lata", 4),
  rep("NE.LibTB", 1),
  rep("NE.NG", 3),
  rep("NE.Youri", 1),
  rep("Tz.PEM", 26),
  rep("Tz.UNG", 21)))

mito_d$pop = factor(c(
  rep("SRA", 3),
  rep("NE", 18),
  rep("TZ_PEM", 26),
  rep("TZ_UNG", 21)))


plot(mito_d$C1, mito_d$C2, col=c(rep("chartreuse3", 3), rep("deeppink", 18), rep("blue", 26), rep("cyan", 21)), lwd=1, pch=19, cex=1.5, xlab="PC 1", ylab="PC 2", main = "MDS of Sh.mito (gatk/plink)")

legend(0.03, 0.025,
       c("SRA", "Niger", "Tanzania (PEM)", "Tanzania (UNG)"),
       col = c("chartreuse3", "deeppink", "blue", "cyan"),
       cex = 0.8,
       lwd = 7, lty = 1)
