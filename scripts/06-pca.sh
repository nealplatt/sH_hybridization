#...............................................................................
#...............................................................................
#PCA - needs LD pruned snps and to be in the plink ped format

#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

mkdir $RESULTS_DIR/pca
cd $RESULTS_DIR/pca

#the goal is to end up with 3 runs of PCA
#1) all samples
#2) curs, bovis, haem
#3) haem only

################################################################################
#1) ALL SAMPLES
plink \
    --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel_LD-25-5-2.vcf \
    --pca \
    --allow-extra-chr \
    --out cohort_snps_schMan_autosomal_panel_LD-25-5-2_pca

#redo the pca but remove all non-haem/bovis/curassoni samples
grep "#" ../build_snp_panel/cohort_snps_schMan_autosomal_panel_LD-25-5-2.vcf  \
    | tail -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >cohort_snps_schMan_autosomal_panel_LD-25-5-2_pca.samples

################################################################################
#2) HAEM_GROUP (includes bovis)

vcftools \
    --vcf  ../build_snp_panel/cohort_snps_schMan_autosomal_panel_LD-25-5-2.vcf \
    --remove-indv ERR103051 \
    --remove-indv ERR119613 \
    --remove-indv ERR310937 \
    --remove-indv ERR310940 \
    --remove-indv ERR539855 \
    --remove-indv ERR539857 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >haem_group_snps_schMan_autosomal_panel_LD-25-5-2_pca.vcf

plink \
    --vcf haem_group_snps_schMan_autosomal_panel_LD-25-5-2_pca.vcf \
    --pca \
    --allow-extra-chr \
    --out haem_group_snps_schMan_autosomal_panel_LD-25-5-2_pca


grep "#" haem_group_snps_schMan_autosomal_panel_LD-25-5-2_pca.vcf \
    | tail -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >haem_group_snps_schMan_autosomal_panel_LD-25-5-2_pca.samples

################################################################################
#3) HAEM_ONLY (includes haem from SRA)
vcftools \
    --vcf haem_group_LD-150-5-02.recode.vcf \
    --remove-indv ERR103048 \
    --recode \
    --recode-INFO-all \
    --out haem_LD-150-5-02

plink \
    --vcf haem_LD-150-5-02.recode.vcf \
    --pca \
    --allow-extra-chr \
    --out haem_LD-150-5-02_pca


grep "#" haem_LD-150-5-02.recode.vcf  \
    | tail -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >haem_LD-150-5-02_pca.samples

################################################################################
# PLOTTING in R

library(scales)

setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/pca")

all_evec<-read.table("all_LD-150-5-02_pca.eigenvec", header=FALSE)
all_eval<-read.table("all_LD-150-5-02_pca.eigenval", header=FALSE)
all_samples<-read.table("all_LD-150-5-02_pca.samples", header=FALSE)

group_evec<-read.table("haem_group_LD-150-5-02_pca.eigenvec", header=FALSE)
group_eval<-read.table("haem_group_LD-150-5-02_pca.eigenval", header=FALSE)
group_samples<-read.table("haem_group_LD-150-5-02_pca.samples", header=FALSE)

haem_evec<-read.table("haem_LD-150-5-02_pca.eigenvec", header=FALSE)
haem_eval<-read.table("haem_LD-150-5-02_pca.eigenval", header=FALSE)
haem_samples<-read.table("haem_LD-150-5-02_pca.samples", header=FALSE)

#haem_sra=orange
#bovis=green
#matt=purple
#inter=yellow
#cur=grey
#marg=brown
#niger=deeppink
#tz_pem=blue
#tz_ung=aquamarine

all_colors<-c(rep("orange", 2), 
          rep("green", 1), 
          rep("purple", 1),
          rep("yellow", 1), 
          rep("grey", 1),
          rep("brown", 1),
          rep("purple", 2),
          rep("orange", 1),
          rep("deeppink", 48), 
          rep("blue", 26),
          rep("aquamarine", 21))

svg("all_LD-150-5-02_pca.svg", width=7, height=7)
plot(all_evec$V3, all_evec$V4, 
     col=alpha(all_colors, 0.4), 
     cex=1.5, 
     pch=19, 
     xlab="PC1 (4.59%)", 
     ylab="PC2 (3.7%)", 
     main="PCA unlinked autosomes"
     )


legend(0.25, 0.25,
      c("S. bovis", 
         "S. curassoni",
         "S. intercalatum", 
         "S. margrebowiei", 
         "S. mattheei", 
         "S. haematobium. (Niger)", 
         "S. h. (Tanzania, Pemba Island)", 
         "S. h. (Tanzania, Ungua Island)",
         "S. h. (Egypt)"), 
       col=alpha(c("green", "grey", "yellow", "brown", "purple", "deeppink", "blue", "aquamarine", "orange"), 0.4),
       cex = 1,
       lwd = 0, 
       lty = 1,
       pch=19
)
dev.off()

svg("all_LD-150-5-02_pca_eigenval.svg", width=7, height=7)
barplot(all_eval$V1, 
        ylim=c(0, 5), 
        xlim=c(0,20), 
        xlab="Principal Component", 
        ylab="Eigenvalue percent", 
        col=c(rep("grey", 2), 
              rep("white", 18))
        
)

dev.off()


