library(scales)
setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/introgression_tracts/")

blocks<-read.table("blocks.tsv", header=T, sep='\t')
dev.off()


#filter
#remove blocks with less than 500 bp
blocks<-blocks[blocks$length_bp>500,]

#remove blocks with less than 500 bp
blocks<-blocks[blocks$snps_in_block>2,]

#remove bovis samples
blocks<-blocks[blocks$indiv_hap!="ERR103048_A",]
blocks<-blocks[blocks$indiv_hap!="ERR103048_B",]

#remove outlier blocks
blocks<-blocks[blocks$length_bp < 3.5e6,]

morgans=(blocks$length_bp/284000)/100

hist(blocks$length_bp, breaks=1500, xlim=c(0,5e6), col="grey", ylim=c(0,2500), main = "Bovis introgression blocks", xlab = "Tract length")
x<-hist(blocks$length_bp)

fit <- nls( x$density ~ T * exp( 1/(T * 3.4e-6)), start=list(T=1));


plot(x=blocks$length_bp, y=blocks$snps_in_block)

max(blocks$snps_in_block)
max(blocks$length_bp)

######################################################################################
ld<-read.table("ne_autapomoprhic.ld", header=T, sep='\t')
distance=(ld$POS2-ld$POS1+1)/287000
plot(x=distance, y=ld$R.2)
plot(x=distance, y=d)
r2<-ld$R.2
D<-ld$D
d<-abs(ld$D)
d
morgan=distance/100
plot(x=morgan, y=d, col="black", pch=19)

fit <- nls( d ~ a * exp( -T * morgan), start=list(a=1,T=1));

summary(fit)
lines( morgan, predict(fit, list(x=morgan)), col='red');
1.7e1

distance

1/((36.88/287000/100) * 0.95) 

mean(distance)

distance
mean(distance)
mean(blocks$length_bp)

lm<-mean(blocks$length_bp/284000/100)
lm

1/(lm*0.95)
lm
