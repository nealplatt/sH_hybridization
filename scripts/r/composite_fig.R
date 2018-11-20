setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/composite")

#load necessary packages
library(scales)
library(binr)

source("https://bioconductor.org/biocLite.R")
biocLite("Sushi")
biocLite("rtracklayer")
library(GenomicRanges)
library(rtracklayer)
library(Sushi)

###########################################################################################
#set up a color palette

anc_colors <- colorRampPalette(c("blue", "gray", "deeppink"))

###########################################################################################
#PCADMIX RESULTS
#read in data
vit<-read.table("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/pcadmix/control-bov-admixed_maf00.vit.txt", sep=" ", header=FALSE, row.names=1)

#add column names to date
colnames(vit)<-c(seq(0, length(vit)-1))


#generate ancestry by averaging the number of "1" haplotypes (bovis) per window
anc<-round(colMeans(vit[0:5484]), digits=2)

#get the snp coordinate (from python script)
coords<-read.table("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/pcadmix/control-bov-admixed_maf00.snp_coords.txt", sep="\t", header=FALSE, row.names=1)

pca_color<-anc_colors(100)[as.numeric(cut(anc,breaks=100))]

pca_color=anc
pca_color[coords$V2=="SM_V7_1"]="black"
pca_color[coords$V2=="SM_V7_2"]="grey"
pca_color[coords$V2=="SM_V7_3"]="black"
pca_color[coords$V2=="SM_V7_4"]="grey"
pca_color[coords$V2=="SM_V7_5"]="black"
pca_color[coords$V2=="SM_V7_6"]="grey"
pca_color[coords$V2=="SM_V7_7"]="black"

#build a color vector to id SNPs with high percentages of bovis ancestry (sig=0.995)
sig<-quantile(anc, 0.995, na.rm = TRUE)

#identify significant snps
pca_color[anc>sig]="red"

###########################################################################################

#plot sprime results
sprime<-read.table("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/sprime/niger_tz-out_auto_sprime_2018-08-13.cul_pos",header = FALSE, sep = "\t")
colnames(sprime)<-c("chr", "site", "cul_site", "score")

sprime_color<-vector()
sprime_color[sprime$chr=="SM_V7_1"]="black"
sprime_color[sprime$chr=="SM_V7_2"]="grey"
sprime_color[sprime$chr=="SM_V7_3"]="black"
sprime_color[sprime$chr=="SM_V7_4"]="grey"
sprime_color[sprime$chr=="SM_V7_5"]="black"
sprime_color[sprime$chr=="SM_V7_6"]="grey"
sprime_color[sprime$chr=="SM_V7_7"]="black"

plot(x=sprime$cul_site, y=sprime$score, 
     ylim=c(0, 250000),
     pch=19,
     col=sprime_color,
     xaxt="n",
     xlab="")


###########################################################################################
#plot fst results
fst<-read.table("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/fst/ne_vs_tz_fst_w500k_50k.tsv", header=TRUE, sep="\t")
colnames(fst)<-c("chr", "start", "stop", "n_var", "weighted_fst", "mean_fst")


chr1_len<-88881357
chr2_len<-48130368
chr3_len<-50458499
chr4_len<-47279781
chr5_len<-25256119
chr6_len<-24989083
chr7_len<-19288021

chr1_cul_len<-0
chr2_cul_len<-chr1_cul_len + chr1_len
chr3_cul_len<-chr2_cul_len + chr2_len
chr4_cul_len<-chr3_cul_len + chr3_len
chr5_cul_len<-chr4_cul_len + chr4_len
chr6_cul_len<-chr5_cul_len + chr5_len
chr7_cul_len<-chr6_cul_len + chr6_len

fst_cul_pos<-vector()
for(i in seq(1, length(fst$start))){
  if ( fst$chr[i] == "SM_V7_1" ) {
    fst_cul_pos[i]<-fst$start[i] +  chr1_cul_len
  } else if ( fst$chr[i] == "SM_V7_2" ) {
    fst_cul_pos[i]<-fst$start[i] +  chr2_cul_len
  } else if ( fst$chr[i] == "SM_V7_3" ) {
    fst_cul_pos[i]<-fst$start[i] +  chr3_cul_len
  } else if ( fst$chr[i] == "SM_V7_4" ) {
    fst_cul_pos[i]<-fst$start[i] +  chr4_cul_len
  } else if ( fst$chr[i] == "SM_V7_5" ) {
    fst_cul_pos[i]<-fst$start[i] +  chr5_cul_len
  } else if ( fst$chr[i] == "SM_V7_6" ) {
    fst_cul_pos[i]<-fst$start[i] +  chr6_cul_len
  } else {
    fst_cul_pos[i]<-fst$start[i] +  chr7_cul_len
  }
  
}

fst_color=vector()
fst_color[fst$chr=="SM_V7_1"]="black"
fst_color[fst$chr=="SM_V7_2"]="grey"
fst_color[fst$chr=="SM_V7_3"]="black"
fst_color[fst$chr=="SM_V7_4"]="grey"
fst_color[fst$chr=="SM_V7_5"]="black"
fst_color[fst$chr=="SM_V7_6"]="grey"
fst_color[fst$chr=="SM_V7_7"]="black"

sig<-quantile(fst$weighted_fst, 0.999, na.rm = TRUE)
fst_color[fst$weighted_fst>sig]="red"

###########################################################################################
#bayescan fst results
bs<-read.table("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/bayescan/test.fst",header = TRUE, sep = "\t")

sites<-read.table("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/bayescan/sites",header = FALSE, sep = "\t")
colnames(sites)<-c("chr", "bp")

bayes_cul_pos<-vector()
for(i in seq(1, length(sites$chr))){
  if ( sites$chr[i] == "SM_V7_1" ) {
    bayes_cul_pos[i]<-sites$bp[i] +  chr1_cul_len
  } else if ( sites$chr[i] == "SM_V7_2" ) {
    bayes_cul_pos[i]<-sites$bp[i] +  chr2_cul_len
  } else if ( sites$chr[i] == "SM_V7_3" ) {
    bayes_cul_pos[i]<-sites$bp[i] +  chr3_cul_len
  } else if ( sites$chr[i] == "SM_V7_4" ) {
    bayes_cul_pos[i]<-sites$bp[i] +  chr4_cul_len
  } else if ( sites$chr[i] == "SM_V7_5" ) {
    bayes_cul_pos[i]<-sites$bp[i] +  chr5_cul_len
  } else if ( sites$chr[i] == "SM_V7_6" ) {
    bayes_cul_pos[i]<-sites$bp[i] +  chr6_cul_len
  } else {
    bayes_cul_pos[i]<-sites$bp[i] +  chr7_cul_len
  }
  
}

bayes_color<-vector()
bayes_color[sites$chr=="SM_V7_1"]="black"
bayes_color[sites$chr=="SM_V7_2"]="grey"
bayes_color[sites$chr=="SM_V7_3"]="black"
bayes_color[sites$chr=="SM_V7_4"]="grey"
bayes_color[sites$chr=="SM_V7_5"]="black"
bayes_color[sites$chr=="SM_V7_6"]="grey"
bayes_color[sites$chr=="SM_V7_7"]="black"

bayes_color[bs$qval<=0.05]="red"

#######################################################################
#xpehh

xpehh<-read.csv("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/rehh/xpehh_maf00_2018-09-06.csv",header = TRUE, sep = ",")
colnames(xpehh)<-c("snp", "chr", "pos", "xpehh", "sig")

chr4_xpehh<-xpehh[xpehh$chr == 4, ]

xpehh_colors=rep("grey", length(chr4_xpehh$xpehh))

xpehh_colors[abs(chr4_xpehh$sig)>=2.301]="red"
#######################################################################
#pi 
pi<-read.csv("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/div/window.pi.csv",header = TRUE, sep = ",")

#######################################################################
#SUSHI
bed<-read.table(file="Sm_v7.0.bed",sep="\t", header = FALSE)

colnames(bed)<-c("chrom", "start", "end", "name", "score", "orient", "prog", "type", "unk", "parent")

genes<-bed[bed$type=="gene", ]


#######################################################################

#define target region on chr4
target_region_start=19.8e6
target_region_end=20.4e6


dev.off()
svg("composite_2018-09-07.svg",width=17.8/2.5, height =12/2.5)

  par(mfrow=c(6,1))
  par(mar=c(0.5,5,0,0))
  par(mgp=c(3,1,0))
  par(oma = c(4, 0.5, 3, 0.5))
  par(las=1)
  
  #----------------------------------------------------------------------------------------------------
  #plot fst
  plot(fst_cul_pos, fst$weighted_fst, 
       pch=19,
       xlim=c(0,305000000),
       ylim=c(0,1),
       col=fst_color,
       ylab=expression(italic("F"["ST"])),
       xaxt="n",
       xlab="",
       xaxs="i",
       yaxp  = c(0, 1, 2),
       cex.lab=1.1,
       cex=0.75
  )
   
  #set axis labels and positions  
  x_axis_labels<-c("Chr 1", "Chr 2", "Chr 3", "Chr 4", "Chr 5", "Chr 6", "Chr 7")
  x_axis_label_pos<-c(chr1_cul_len+(chr1_len/2), 
                      chr2_cul_len+(chr2_len/2), 
                      chr3_cul_len+(chr3_len/2), 
                      chr4_cul_len+(chr4_len/2), 
                      chr5_cul_len+(chr5_len/2), 
                      chr6_cul_len+(chr6_len/2), 
                      chr7_cul_len+(chr7_len/2)
  )
  #add the tick marks
  axis(side=3, 
       at = c(0, 88881357, 137011725, 187470224, 234750005, 260006124, 284995207, 304283228),
       labels=c(rep("", 8)))
  
  #and labels
  mtext(x_axis_labels, 
        side=3,
        at=x_axis_label_pos,
        line=1,
        cex=0.75
  )
  
  #-BAYESCAN-----------------------------
  plot(x=cul_pos, y=bs$alpha, 
       pch=19,
       ylim=c(-0.1,2),
       xlim=c(0,305000000),
       ylab=expression(alpha), 
       col=bayes_color,
       xaxt="n",
       xaxs="i",
       cex.lab=1.3,
       cex = 0.75,
       yaxp  = c(0, 2, 2)
  )
  
  
  #-ANC-----------------------------
  #back to standard margin
  #par(mar=c(0.5,5,0,0))
  
  plot(x=coords$V3, y=anc, 
       pch=19, 
       xlim=c(0,305000000),
       ylim=c(0,1),
       col=pca_color,
       ylab="Ancestry",
       xaxt="n",
       xlab="",
       xaxs="i",
       xaxt="n",
       cex.lab=1.1,
       cex=0.75,
       yaxp  = c(0, 1, 2)
      )
  
  #-PI--------------------------------
  #ZOOMED IN
  plot(x = pi$MID_POS, 
       y = pi$RATIO, 
       pch=19, 
       xlim=c(target_region_start, target_region_end),
       type="l", 
       lwd=2,
       ylab=expression(paste("log(", pi["NE"]/pi["TZ"], ")")),
       xlab="",
       xaxt="n",
       xaxs="i",
       cex.lab=0.9,
       yaxp  = c(-2, 2, 2)
  )
  
  #-XPEHH-----------------------------
  #ZOOMED IN
  plot(x = chr4_xpehh$pos, 
       y = chr4_xpehh$xpehh, 
       pch=19, 
       xlim=c(target_region_start, target_region_end),
       ylim=c(-7, 7),
       type="l", 
       lwd=2,
       ylab="xpEHH",
       xlab="",
       cex.lab=1.1,
       yaxp = c(-6, 6, 2)
  )
  
  #-SUSHI--------------------------------
  #ZOOMED IN
  plotBed(beddata = genes,
          chrom = "SM_V7_4",
          chromstart = target_region_start,
          chromend =target_region_end,
          row=3,
          wiggle=0
  )
  labelgenome("SM_V7_4", target_region_start,target_region_end,n=10,scale="Mb")
  
dev.off()


