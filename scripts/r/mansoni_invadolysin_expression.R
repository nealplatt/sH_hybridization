install.packages("viridis")
library(viridis)

setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/mansoni_expression")
d<-read.csv("mansoni_invadolysin_expression.csv", header=TRUE)

colors<-c("white", "green")


rbPal <- colorRampPalette(c('white','green'))

rbPal(10)[as.numeric(cut(d[1,3:13],breaks = 10))]

image(1:nrow(d), 1:ncol(d), as.matrix(d), col=cols)

dat <- sapply( dat, as.numeric )


aa<-read.csv("invadolysin_AAs_matrix_no_headers.csv", sep = ",", header=FALSE)
aa

cols <- c(
  'F' = "grey",
  'V' = "grey",
  'T' = "grey",
  'K' = "grey",
  'I' = "grey",
  'H' = "grey",
  'Y' = "grey",
  'D' = "grey",
  'G' = "grey",
  '0' = "white"
)

cols <- c(
  '1' = "cornsilk",
  '0' = "grey28"
)


image(1:nrow(aa), 1:ncol(aa), as.matrix(aa), col=cols)










data<-dat[,3:13]
col_grad <- colorRampPalette(c("white", "green4"))(n = 25)
heatmap(dat[,3:13], 
        Colv = NA, 
        Rowv = NA, 
        scale = "row", 
        col=col_grad,  
        labRow = c("Smp_167070.2",
                   "Smp_167070.1",
                   "Smp_336180.1",
                   "Smp_331850.1",
                   "Smp_314000.1",
                   "Smp_303070.1",
                   "Smp_303760.1",
                   "Smp_303760.2",
                   "Smp_247860.1",
                   "Smp_247870.1",
                   "Smp_336180.5",
                   "Smp_336180.3",
                   "Smp_127030.1",
                   "Smp_247850.1",
                   "Smp_334410.1",
                   "Smp_153930.1"),
        labCol = c("Miracidia",
                   "Sporocyst [48h]",
                   "Cercariae",
                   "Schistosomula [3h]",
                   "Schistosomula [24h]",
                   "Adult Male [21d]",
                   "Adult Female [21d]",
                   "Adult Male [28d]",
                   "Adult Female [28d]",
                   "Adult Male [38d]",
                   "Adult Female [38d]"),
        
       )

viridis(10)
library(gplots)
svg("mansoni_expression.svg")
heatmap.2(dat[,3:13],
          trace="none",
          scale="row",
          Rowv = FALSE,
          Colv = FALSE,
          col=viridis(10),
          key.title=NA,
          key.xlab=NA,
          dendrogram="none",
          density.info="none",
          margins=c(12,9),
          labRow = c("Smp_167070.2",
                     "Smp_167070.1",
                     "Smp_336180.1",
                     "Smp_331850.1",
                     "Smp_314000.1",
                     "Smp_303070.1",
                     "Smp_303760.1",
                     "Smp_303760.2",
                     "Smp_247860.1",
                     "Smp_247870.1",
                     "Smp_336180.5",
                     "Smp_336180.3",
                     "Smp_127030.1",
                     "Smp_247850.1",
                     "Smp_334410.1",
                     "Smp_153930.1"),
          labCol = c("Miracidia",
                     "Sporocyst [48h]",
                     "Cercariae",
                     "Schistosomula [3h]",
                     "Schistosomula [24h]",
                     "Adult Male [21d]",
                     "Adult Female [21d]",
                     "Adult Male [28d]",
                     "Adult Female [28d]",
                     "Adult Male [38d]",
                     "Adult Female [38d]"),
          symm=F,
          symkey=F,
          symbreaks=FALSE,
          colsep=c(1,2,3,5),
          keysize = 1
          )
dev.off()

colSums(data)
barplot(log10(colSums(data)), ylim=c(0,3))

rowSums(data)
barplot(log10(rowSums(data)), xlim=c(0,3), horiz=TRUE)

image(1:10, 1, as.matrix(1:10), 
      +       col=viridis(10),
      +       xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")


dev.off()
barplot(data[13,])

x<-read.csv("reduced_expression.csv", header=TRUE)
barplot(as.matrix(x), 
        beside=TRUE, 
        col=c(rep("grey", 10), 
              rep(c("blue", "pink"), 3)
              ), 
        ylab="TPM", 
        ylim=c(0, 7),
        legend=c(rownames(x)) 
        )
        
