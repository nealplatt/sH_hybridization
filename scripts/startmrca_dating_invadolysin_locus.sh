#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# startmrca_dating_invadolysin_locus.sh - use startmrca to date time since
#   common ancestor of the incadolysin locus 

# Uses a conda to manage the enviroment and is run on a high memory node; 125gb
#Set up the environment
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir tmrca

cd tmrca


grep "#" ../beagle/auto_beagle_maf05.vcf | tail -n 1 | cut -f13- | sed 's/\t/\n/g' >samples.list
grep Sh.NE samples.list >ne.list
grep Sh.TZ samples.list >tz.list


#SM_V7_4 20033013        KL250964.1:40108 

vcftools \
    --vcf ../beagle/auto_beagle_maf05.vcf \
    --keep ../samples.list \
    --chr SM_V7_4 \
    --from-bp 19033012 \
    --to-bp 21033014 \
    --recode \
    --stdout \
    | gzip \
        >target_snp.vcf.gz

#get probed regions
${ENVIRONMENTS["TITAN SINGULARITY"]} \
    /usr/software/progressiveCactus/submodules/hal/bin/halLiftover  \
        $RESULTS_DIR/wga/schMan7_vs_schMan1.hal \
        schHae_v1 \
        ../../data/schHae_v1_probes.bed \
        schMan_v7 \
        schHae_v1_probes_schMan_lifted.bed

bedtools sort \
    -i schHae_v1_probes_schMan_lifted.bed \
    | grep "SM_V7_4\s" \
        >SM_V7_4_probes.bed


awk -v OFS='\t' {'print $1,$2'} ../../data/genome/schMan_v7.fa.fai > schMan_v7.genome

cat schMan_v7.genome | grep "SM_V7_4\s" >SM_V7_4.genome


bedtools \
    complement \
    -i SM_V7_4_probes.bed \
    -g SM_V7_4.genome \
    >unprobed_regions.bed

#This portion of the script is run in R
R

library(startmrca)

for (chain in c(1:10)){
  run.startmrca(vcf.file    = "target_snp.vcf.gz",
              rec.file      = NULL, 
              sample.ids    = "ne.list", 
              refsample.ids = "tz.list",
              mut.rate      = 8.1e-9,
              rec.rate      = 3.4e-8, 
              nsel          = 46,
              nanc          = 47,
              chain.length  = 50000,
              nanc.post     = 1000,
              pos           = 20033013,
              sel.allele    = 0,
              proposal.sd   = 20,
              bed.file      ="unprobed_regions.bed")
    
    file_name=paste("mcmc.target_t.chain", chain, ".csv", sep="")
    load("target_snp_mcmc_list.RDATA")
    write.csv(mcmc.output$t.chain, file=file_name)

    #chain<-mcmc.output$t.chain[,1]
    #post_burnin<-chain[1500:2000]
    #samples=post_burnin[1:100*5]
    #quantile(samples, c(0.025, 0.975))
}

# in bash (burnin and thin)
for CHAIN in $(seq 1 10); do
    tail -n 10000 mcmc.target_t.chain$CHAIN.csv \
    | awk 'NR%10==1' \
        >mcmc.target_t.chain$CHAIN"_burnin.csv"
done

