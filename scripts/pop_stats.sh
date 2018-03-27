source /master/nplatt/sH_hybridization/scripts/set-env.sh

cd /master/nplatt/sH_hybridization/results/pop_stats

#make populations lists for each location of interest
mkdir lists
#<make lists>
# $> ls lists/
# Sh.NE.list  Sh.TZ.list  Sh.TZ_PEM.list  Sh.TZ_UNG.list

vcftools \
    --vcf ../filter_cohort_vcf/sHaem_filtered.vcf \
    --weir-fst-pop lists/Sh.NE.list \
    --weir-fst-pop lists/Sh.TZ.list \
    --out NE_vs_TZ



# filter "linked" snps for population assignment/pca
# see file:///C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/docs/congenomics_plink_tutorial_davey.pdf
plink \
    --vcf ../filter_cohort_vcf/sHaem_filtered.vcf \
    --allow-extra-chr \
    --indep-pairwise 25 5 0.2 \
    --out LDprunedout

singularity shell gatk SelectVariants -V sHaem_filtered_cohort_named.vcf -R $REFERENCE -O sHaem_filtered_cohort_named_LDpruned.vcf --exclude-ids LDprunedout.prune.out



plink --vcf sHaem_filtered_cohort_named.vcf --allow-extra-chr --recode structure --out sHaem_filtered_cohort_named



#10 is the number of replicates at each value of K.
#check out param files and update
for rep in {1..10} 
do 
    structure -m mainparams -K 2 -i sHaem_filtered_cohort_named -o outfile_k${SGE_TASK_ID}_rep${rep} <change burnin and reps>
done

#define INFILE sHaem_filtered_cohort.recode.strct_in
#define OUTFILE sHaem_filtered_cohort.structure
#define NUMINDS 89
#define NUMLOCI 133215
#define LABEL 1
#define POPDATA 1
#define POPFLAG 0
#define PHENOTYPE 0
#define EXTRACOLS 0
#define PHASEINFO 0
#define MISSING -9
#define PLOIDY 2
#define ONEROWPERIND 1
#define MARKERNAMES 1 
#define MAPDISTANCES 1
#define MAXPOPS 2
#define BURNIN 100000
#define NUMREPS 1000000

