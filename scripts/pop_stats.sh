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

mv LDprunedout.prune.out LDprunedout.prune.out.list

singularity exec ../../snpCalling_v0.0.7.img \
    gatk SelectVariants \
        --exclude-ids LDprunedout.prune.out.list \
        -V ../filter_cohort_vcf/sHaem_filtered.vcf \
        -R $REFERENCE \
        -O sHaem_filtered_LDpruned.vcf \


vcftools \
    --vcf sHaem_filtered_LDpruned.vcf \
    --weir-fst-pop lists/Sh.NE.list \
    --weir-fst-pop lists/Sh.TZ.list \
    --out NE_vs_TZ-LDpruned

#admixture
plink --vcf sHaem_filtered_LDpruned.vcf --out sHaem_filtered_LDpruned --recode12 --allow-extra-chr

wget https://www.genetics.ucla.edu/software/admixture/binaries/admixture_linux-1.3.0.tar.gz

for K in {1..12}; do
    admixture_linux-1.3.0/admixture --cv=100 sHaem_filtered_LDpruned.ped $K >$K.log &
done

#do CV plot and plot Q values
#split this up by NE and TZ samples

#select TZ
singularity exec ../../snpCalling_v0.0.7.img \
    gatk SelectVariants \
        --se 'Sh.TZ_' \
        -V sHaem_filtered_LDpruned.vcf \
        -R $REFERENCE \
        -O Sh.TZ_filtered_LDpruned.vcf

#select NE
singularity exec ../../snpCalling_v0.0.7.img \
    gatk SelectVariants \
        --se 'Sh.NE_' \
        -V sHaem_filtered_LDpruned.vcf \
        -R $REFERENCE \
        -O Sh.NE_filtered_LDpruned.vcf


plink --vcf Sh.TZ_filtered_LDpruned.vcf --out Sh.TZ_filtered_LDpruned --recode12 --allow-extra-chr
plink --vcf Sh.NE_filtered_LDpruned.vcf --out Sh.NE_filtered_LDpruned --recode12 --allow-extra-chr

for K in {1..12}; do
    admixture_linux-1.3.0/admixture --cv=100 Sh.NE_filtered_LDpruned.ped $K >NE.$K.log &
done

wait

for K in {1..12}; do
    admixture_linux-1.3.0/admixture --cv=100 Sh.TZ_filtered_LDpruned.ped $K >TZ.$K.log &
done


#PCA
plink --file sHaem_filtered_LDpruned --genome --allow-no-sex --allow-extra-chr -out sHaem_filtered_LDpruned

plink --file sHaem_filtered_LDpruned \
    --read-genome sHaem_filtered_LDpruned.genome \
    --cluster --mds-plot 2 \
    --allow-extra-chr

#R code
d <- read.table("plink.mds", h=T)

d$specLoc = factor(c(rep("NE.Dai", 13),
                 rep("NE.Doki", 1),
                 rep("NE.Kar", 8),
                 rep("NE.Lata", 6),
                 rep("NE.LibTB", 5),
                 rep("NE.NG", 4),
                 rep("NE.Seb", 3),
                 rep("NE.Tiag", 1),
                 rep("NE.YK", 4),
                 rep("NE.Youri", 2),
                 rep("Tz.PEM", 25),
                 rep("Tz.UNG", 19)))

d$pop = factor(c(rep("NE", 47),
                 rep("TZ", 44)))


plot(d$C1, d$C2, col=c(rep("deeppink", 47), rep("blue", 44)), pch=19, xlab="PC 1", ylab="PC 2", main = "MDS of Sh (gatk/plink)")

#legend(x=0.1, y=0.3, c("NE.Dai", "NE.Doki", "NE.Kar", "NE.Lata", "NE.LibTB", "NE.NG", "NE.Seb", "NE.Tiag", "NE.Yk", "NE.Youri", "TZ.PEM", "TZ.UNG"), pch=19, col=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))

legend(0.03, 0.025,
       c("Niger", "Tanzania"),
       col = c("deeppink", "blue"),
       cex = 0.8,
       lwd = 1, lty = 1)


