setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/f3")
af3_ne_tz_bov=scan(file="af3_ne_tz_bov")
af3_ne_tz_curs=scan(file="af3_ne_tz_curs")
f3_ne_tz_bov=scan(file="f3_ne_tz_bov")
f3_ne_tz_curs=scan(file="f3_ne_tz_curs")

af3_tz_ne_bov=scan(file="af3_tz_ne_bov")
af3_tz_ne_curs=scan(file="af3_tz_ne_curs")
f3_tz_ne_bov=scan(file="f3_tz_ne_bov")
f3_tz_ne_curs=scan(file="f3_tz_ne_curs")


#bovis=green
#cur=grey


svg("f3_2018-09-14.svg")
par(mar=c(5,10,5,5))
boxplot(af3_ne_tz_bov, af3_ne_tz_curs, f3_ne_tz_bov, f3_ne_tz_curs, f3_tz_ne_bov, f3_tz_ne_curs, 
        names =c("aF3(NE: TZ, BOV)", "aF3(NE: TZ, CURS)", "F3(NE: TZ, Bov)", "F3(NE: TZ, CURS)", "F3(TZ: NE, Bov)", "F3(TZ: NE, Curs)"),
        col="light grey",
        pch=19, 
        cex=0.5, 
        lty=1,
        xlab="f3",
        horizontal = TRUE,
        notch=TRUE,
        outline=FALSE,
        boxwex=0.3,
        staplewex=0,
        las=2,
        lwd=1.5
        
      
)
abline(v = 0, col = "red", lty = 3, lwd=2) 
dev.off()

boxplot(af3_ne_tz_bov, af3_ne_tz_curs)
