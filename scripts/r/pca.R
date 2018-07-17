library(scales)

setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/pca")

#color codes
#haem_sra=orange
#bovis=green
#matt=purple
#inter=yellow
#guin=coral2
#cur=grey
#marg=brown
#niger=deeppink
#tz_pem=blue
#tz_ung=aquamarine

#------------------------------------------------------------------
#plot haem only samples
haem_evec<-read.table("haemOnly_auto_maf_ld_pca.eigenvec", header=FALSE)
haem_eval<-read.table("haemOnly_auto_maf_ld_pca.eigenval", header=FALSE)
haem_samples<-read.table("haemOnly_auto_maf_ld_pca.samples", header=FALSE)

haem_colors<-c(
  rep("orange", 1),
  rep("deeppink", 48), 
  rep("blue", 26),
  rep("aquamarine", 21))

svg("haemOnly_auto_maf_ld_pca_2018-07-06.svg", width=7, height=7)
plot(haem_evec$V3, haem_evec$V4, 
     col=alpha(haem_colors, 0.7), 
     cex=1.5, 
     pch=19, 
     xlab="PC1 (22.1%)", 
     ylab="PC2 (2%)", 
     main="PCA unlinked autosomes"
)

legend(-.1, 0.8,
       c("S. haematobium. (Niger)", 
         "S. h. (Tanzania, Pemba Island)", 
         "S. h. (Tanzania, Ungua Island)",
         "S. h. (Egypt)"), 
       col=alpha(c("deeppink", "blue", "aquamarine", "orange"), 0.7),
       cex = 1,
       lwd = 0, 
       lty = 1,
       pch=19
)
dev.off()
#------------------------------------------------------------------


#------------------------------------------------------------------
#plot haem and bovis samples
group_evec<-read.table("haem_auto_maf_ld_pca.eigenvec", header=FALSE)
group_eval<-read.table("haem_auto_maf_ld_pca.eigenval", header=FALSE)
group_samples<-read.table("haem_auto_maf_ld_pca.samples", header=FALSE)

group_colors<-c(
  rep("orange", 1),
  rep("green", 1), 
  rep("deeppink", 48), 
  rep("blue", 26),
  rep("aquamarine", 21))

svg("haem_auto_maf_ld_pca_2018-07-06.svg", width=7, height=7)
plot(group_evec$V3, group_evec$V4, 
     col=alpha(group_colors, 0.7), 
     cex=1.5, 
     pch=19, 
     xlab="PC1 (22.5%)", 
     ylab="PC2 (8.9%)", 
     main="PCA unlinked autosomes"
)

legend(-0.1, -0.6,
       c("S. bovis", 
         "S. haematobium. (Niger)", 
         "S. h. (Tanzania, Pemba Island)", 
         "S. h. (Tanzania, Ungua Island)",
         "S. h. (Egypt)"), 
       col=alpha(c("green", "deeppink", "blue", "aquamarine", "orange"), 0.7),
       cex = 1,
       lwd = 0, 
       lty = 1,
       pch=19
)
dev.off()
#------------------------------------------------------------------



#------------------------------------------------------------------
#plot all samples

all_evec<-read.table("auto_maf_ld_pca.eigenvec", header=FALSE)
all_eval<-read.table("auto_maf_ld_pca.eigenval", header=FALSE)
all_samples<-read.table("auto_maf_ld.samples", header=FALSE)

all_colors<-c(rep("orange", 1), 
              rep("green", 1), 
              rep("purple", 1),
              rep("yellow", 1),
              rep("coral2", 1),
              rep("grey", 1),
              rep("brown", 1),
              rep("purple", 2),
              rep("orange", 1),
              rep("deeppink", 48), 
              rep("blue", 26),
              rep("aquamarine", 21))


svg("auto_maf_ld_pca_2018-07-06.svg", width=7, height=7)
plot(all_evec$V3, all_evec$V4, 
     col=alpha(all_colors, 0.7), 
     cex=1.5, 
     pch=19, 
     xlab="PC1 (30.5%)", 
     ylab="PC2 (18.7%)", 
     main="PCA unlinked autosomes"
     )


legend(0.1, 0.1,
      c("S. bovis", 
         "S. curassoni",
         "S. guiniensis",
         "S. intercalatum", 
         "S. margrebowiei", 
         "S. mattheei", 
         "S. haematobium. (Niger)", 
         "S. h. (Tanzania, Pemba Island)", 
         "S. h. (Tanzania, Ungua Island)",
         "S. h. (Egypt)"), 
       col=alpha(c("green", "grey","coral2",  "yellow", "brown", "purple", "deeppink", "blue", "aquamarine", "orange"), 0.7),
       cex = 1,
       lwd = 0, 
       lty = 1,
       pch=19
)
dev.off()
#------------------------------------------------------------------

svg("all_LD-150-5-02_pca_eigenval.svg", width=7, height=7)
barplot(all_eval$V1, 
        ylim=c(0, 35), 
        xlim=c(0,20), 
        xlab="Principal Component", 
        ylab="Eigenvalue percent", 
        col=c(rep("grey", 2), 
              rep("white", 18))
        
)

dev.off()
