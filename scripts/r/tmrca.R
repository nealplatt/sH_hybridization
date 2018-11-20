library(startmrca)

setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/tmrca")

run.startmrca(vcf.file      = "target_snp.vcf",
              rec.file      = NULL, 
              sample.ids    = "ne.list", 
              refsample.ids = "tz.list",
              mut.rate      = 8.1e-9,
              rec.rate      = 3.4e-8, 
              nsel          = 50,
              nanc          = 20,
              chain.length  = 20,
              nanc.post     = 10,
              pos           = 20033013,
              sel.allele    = 0,
              bed.file      ="unprobed_regions.bed")

load("target_snp_mcmc_list.RDATA")

mcmc.output$anchap.post

chain<-mcmc.output$t.chain[,1]
post_burnin<-chain[1500:2000]
samples=post_burnin[1:100*5]
quantile(samples, c(0.025, 0.975))

for (chain in c(1:10)){
  print(chain)
}


data<-read.csv("mcmc.target_t.chain1.csv")
chain<-data[,2]
post_burnin<-chain[1500:2000]
samples=post_burnin[1:100*5]
quantile(samples, c(0.025, 0.975))
