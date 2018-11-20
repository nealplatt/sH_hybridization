setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/bayescan")
source("plot_R.r")


plot_bayescan("tz-ne_maf05_bi_bayescan_pr10_pi10k_bin50K_ngen50_npb50_thin20_chain001_fst.txt",FDR=0.05)

bs<-read.table("tz-ne_maf05_bi_bayescan_pr10_pi10k_bin50K_ngen50_npb50_thin20_chain001_fst.csv",header = TRUE, sep = ",")

plot(bs$alpha, pch=19, ylim=c(-2,2), ylab="Alpha", xlab="Chromosome")



library(coda)
chain <- read.table(file="test.sel",header=TRUE)
chain<-chain[-c(1)]

chain<-mcmc(chain,thin=10)
plot(chain)
summary(chain)

chain<-mcmc(chain,thin=1000)
summary(chain)


geweke.diag(chain, frac1=0.1, frac2=0.5)

autocorr.diag(chain)


chain <- mcmc(chain,thin=200)
summary(chain)
autocorr.diag(chain)
effectiveSize(chain)
geweke.diag(chain, frac1=0.1, frac2=0.5)


heidel.diag(chain, eps=0.1, pvalue=0.05)

geweke.plot(chain, frac1=0.1, frac2=0.5)
plot(chain)

plot_bayescan("test_pr50_fst.txt",FDR=0.01)
