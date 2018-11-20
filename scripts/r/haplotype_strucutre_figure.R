setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/distribution")

random<-read.table("all.bed", header=FALSE)
actual<-read.table("merged_sizes.bed", header=FALSE)

colnames(random)<-c("chr", "start", "stop", "size")
colnames(actual)<-c("chr", "start", "stop", "size")

t.test(random$size, actual$size)

library(MASS)

poisson.test(c(n1, n2), c(t1, t2), alternative = c("two.sided"))

hist(actual$size, breaks=15)
hist(random$size, breaks=15)

max(actual$size)
max(random$size)


setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/invadolysin_phylo")
geno<-read.csv("invadolysin_genotypes_matrix_no_headers.csv", sep = ",")

cols <- c(
  '0' = "blue",
  '1' = "red",
  '2' = "white"
)
?image()

image(1:nrow(geno), 1:ncol(geno), as.matrix(geno), col=cols)

d<-read.table(text="
0  1  0  3
3  2  1  4
4  1  0  2
3  3  0  1")

cols <- c(
  '0' = "white",
  '0.5' = "grey",
  '2' = "black",
  '3' = "#33CC00",
  '4' = "#009900"
)

image(1:nrow(d), 1:ncol(d), as.matrix(d), col=cols)


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

