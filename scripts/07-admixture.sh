#...............................................................................
#...............................................................................
#...............................................................................
#ADMIXTURE - needs LD pruned snps and to be in the plink ped format

#NEED TO RE-RUN AND EXCLUDE SEX SNPs

#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

mkdir $RESULTS_DIR/admixture
cd $RESULTS_DIR/admixture

#using the LD filtered ped files from the PCA analyses.  goal is to run 2
# analyese.
#1) with all samples
#2) with haematobium, bovis, and curassoni samples.

#consider doing supervised an unsupervised

#ALL SAMPLES
plink \
    --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel_LD-25-5-2.vcf \
    --out cohort_snps_schMan_autosomal_panel_LD-25-5-2 \
    --recode12 \
    --allow-extra-chr

grep "#" ../build_snp_panel/cohort_snps_schMan_autosomal_panel_LD-25-5-2.vcf \
    | head -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >cohort_samples.list

#UNSUPERVISED 
#submit jobs to the cluster
for K in $(seq 1 20); do
    CMD="$WORK_DIR/scripts/admixture_linux-1.3.0/admixture \
        --cv=1000 \
        -j12 \
        cohort_snps_schMan_autosomal_panel_LD-25-5-2.ped \
        $K \
        >all_cv_k$K.log"

    echo $CMD | $QSUB -N all_cv_k$K -o all_cv_k$K.stdout -pe mpi 12
done

#when all runs are done get CV scores
grep CV all_cv_k*.log \
    | sort -n \
    | awk '{print $3"\t"$4}' \
    | sed 's/(//' \
    | sed 's/)//' \
    | sed 's/://' \
    >all_cv_table.tsv


################################################################################
#2) HAEM_GROUP (includes bovis)

#remove all samples not in haem group
vcftools \
    --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel_LD-25-5-2.vcf \
    --remove-indv ERR103051 \
    --remove-indv ERR310940 \
    --remove-indv ERR539855 \
    --remove-indv ERR539857 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >haem_group_snps_schMan_autosomal_panel_LD-25-5-2.vcf

grep "#" haem_group_snps_schMan_autosomal_panel_LD-25-5-2.vcf \
    | head -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >haem_group_samples.list

#convert to ped
plink \
    --vcf haem_group_snps_schMan_autosomal_panel_LD-25-5-2.vcf \
    --out haem_group_snps_schMan_autosomal_panel_LD-25-5-2 \
    --recode12 \
    --allow-extra-chr

#submit jobs to the cluster
for K in $(seq 1 20); do
    CMD="$WORK_DIR/scripts/admixture_linux-1.3.0/admixture \
        --cv=1000 \
        -j12 \
        haem_group_snps_schMan_autosomal_panel_LD-25-5-2.ped \
        $K \
        >group_cv_k$K.log"

    echo $CMD | $QSUB -N group_cv_k$K -o group_cv_k$K.stdout -pe mpi 12
done

#get CV scores when jobs are completed 
grep CV group_cv_k*.log \
    | sort -n \
    | awk '{print $3"\t"$4}' \
    | sed 's/(//' \
    | sed 's/)//' \
    | sed 's/://' \
    >group_cv_table.tsv


#-------------------------------------------------------------------------------
#PLOTTING IN R

library("RColorBrewer")    

setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/admixture")
Q2=read.table("cohort_snps_schMan_autosomal_panel_LD-25-5-2.2.Q")
Q3=read.table("cohort_snps_schMan_autosomal_panel_LD-25-5-2.3.Q")
Q4=read.table("cohort_snps_schMan_autosomal_panel_LD-25-5-2.4.Q")
Q5=read.table("cohort_snps_schMan_autosomal_panel_LD-25-5-2.5.Q")
Q6=read.table("cohort_snps_schMan_autosomal_panel_LD-25-5-2.6.Q")
Q7=read.table("cohort_snps_schMan_autosomal_panel_LD-25-5-2.7.Q")
Q8=read.table("cohort_snps_schMan_autosomal_panel_LD-25-5-2.8.Q")

cohort_names<-scan("cohort_samples.list", what = "character" )

bovis<-"green"
haem_sra<-"orange"
haem_tz_ung<-"blue"
haem_niger<-"deeppink"
haem_tz_pem<-"aquamarine"
matt<-"purple"
marg<-"brown"
inter<-"yellow"
cur<-"grey"

svg("cohort_snps_schMan_autosomal_panel_LD-25-5-2_admixture.svg")
par(mar=c(1,1,1,1))
par(mfrow=c(4,1))
par(cex = 0.001)
barplot(t(as.matrix(Q2)), col=c(bovis, haem_niger), ylab="Ancestry", space=0, border=NA)
barplot(t(as.matrix(Q3)), col=c(haem_tz_ung, haem_niger, bovis), ylab="Ancestry", space=0, border=NA)
barplot(t(as.matrix(Q4)), col=c(matt, haem_niger, bovis, haem_tz_ung), ylab="Ancestry", space=0, border=NA)
barplot(t(as.matrix(Q5)), col=c(haem_tz_pem, haem_niger, bovis, haem_tz_ung, haem_sra), ylab="Ancestry", space=0, border=NA)
barplot(t(as.matrix(Q6)), col=c(haem_tz_pem, bovis, marg, haem_tz_ung, haem_sra, haem_niger), ylab="Ancestry", space=0, border=NA)
dev.off()




# SUPERVISED?


#supervised admixture
#make pop file with haem, curs, bov
for K in $(seq 1 10); do
    CMD="$WORK_DIR/scripts/admixture_linux-1.3.0/admixture \
        --j12 \
        --cv=1000 \
        $GENO_ALL_DIR/all_schisto_smancoord_filtered_LD-25-5-2.ped \
        --supervised
        $K \
        >$GENO_ALL_DIR/cv_k$K.log"

    echo $CMD | $QSUB -N admixture_k$K -o cv_k$K.stdout -pe mpi 12
done

##get CV scores
#grep CV $GENO_ALL_DIR/cv_k*.log \
#    | sort -n \
#    | awk '{print $3"\t"$4}' \
#    | sed 's/(//' \
#    | sed 's/)//' \
#    | sed 's/://' \
#    >$GENO_ALL_DIR/admixture_cv_table.tsv
   

#plot in R on local computer using the script <R SCRIPT>
