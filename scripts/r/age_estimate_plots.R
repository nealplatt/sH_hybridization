setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/introgression_tracts/")

ld_data<-read.table("age_estimates.txt", header = TRUE)

hist(ld_data$est_age_in_gens, xlim=c(0,625), breaks=20, col="red")

setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/tmrca/")

chain1_data<-read.csv("mcmc.target_t.chain1_burnin.csv", header = FALSE)
chain2_data<-read.csv("mcmc.target_t.chain2_burnin.csv", header = FALSE)
chain3_data<-read.csv("mcmc.target_t.chain3_burnin.csv", header = FALSE)
chain4_data<-read.csv("mcmc.target_t.chain4_burnin.csv", header = FALSE)
chain5_data<-read.csv("mcmc.target_t.chain5_burnin.csv", header = FALSE)
chain6_data<-read.csv("mcmc.target_t.chain6_burnin.csv", header = FALSE)
chain7_data<-read.csv("mcmc.target_t.chain7_burnin.csv", header = FALSE)
chain8_data<-read.csv("mcmc.target_t.chain8_burnin.csv", header = FALSE)
chain9_data<-read.csv("mcmc.target_t.chain9_burnin.csv", header = FALSE)
chain10_data<-read.csv("mcmc.target_t.chain10_burnin.csv", header = FALSE)

chain1_data

boxplot(chain1_data$V2,
        chain2_data$V2,
        chain3_data$V2,
        chain4_data$V2,
        chain5_data$V2,
        chain6_data$V2,
        chain7_data$V2,
        chain8_data$V2,
        chain9_data$V2,
        chain10_data$V2,
        ld_data$est_age_in_gens,
        ylim=c(0,700))


vioplot(chain1_data$V2,
        chain2_data$V2,
        chain3_data$V2,
        chain4_data$V2,
        chain5_data$V2,
        chain6_data$V2,
        chain7_data$V2,
        chain8_data$V2,
        chain9_data$V2,
        chain10_data$V2,
        ld_data$est_age_in_gens,
        ylim=c(0,700))

library(vioplot)


blocks<-read.table("blocks_filtered..txt", header = TRUE)
sorted_blocks<-sort(blocks$length_bp, decreasing = TRUE)
plot(sorted_blocks, pch=19)
     