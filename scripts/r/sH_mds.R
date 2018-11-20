d <- read.table("plink.mds", h=T)

d$specLoc = factor(c(rep("NE.Dai", 13),
                 rep("NE.Doki", 1),
                 rep("NE.Kar", 8),
                 rep("NE.Lata", 6),
                 rep("NE.LibTB", 5),
                 rep("NE.NG", 4),
                 rep("NE.Seb", 3),
                 rep("NE.Tiag", 1),
                 rep("NE.YK", 4),
                 rep("NE.Youri", 2),
                 rep("Tz.PEM", 25),
                 rep("Tz.UNG", 19)))

d$pop = factor(c(rep("NE", 47),
                 rep("TZ_PEM", 25),
                 rep("TZ_UNG", 19)))


plot(d$C1, d$C2, 
     col=c(rep("deeppink", 47), rep("blue", 25), rep("green", 19)), 
     #pch=c(rep(1, 47), rep(2, 25), rep(3, 19)), 
     xlab="PC 1", 
     ylab="PC 2", 
     main = "PCA Sh")

legend(x=0.1, y=0.3, c("NE.Dai", "NE.Doki", "NE.Kar", "NE.Lata", "NE.LibTB", "NE.NG", "NE.Seb", "NE.Tiag", "NE.Yk", "NE.Youri", "TZ.PEM", "TZ.UNG"), pch=19, col=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))

legend(0.03, 0.025,
       c("Niger", "Tanzania (Unguja)", "Tanzania (Pemba)"),
       col = c("deeppink", "blue", "green"),
       cex = 0.8,
       lwd = 1, lty = 1)
