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
    --vcf ../build_snp_panel/auto_maf_ld.vcf \
    --out auto_maf_ld \
    --recode12 \
    --allow-extra-chr

grep "#" ../build_snp_panel/auto_maf_ld.vcf \
    | tail -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >samples.list

#UNSUPERVISED 
#submit jobs to the cluster
for K in $(seq 1 20); do
    CMD="$WORK_DIR/scripts/admixture_linux-1.3.0/admixture \
        --cv=1000 \
        -j12 \
        auto_maf_ld.ped \
        $K \
        >all_cv_maf05_k$K.log"

    echo $CMD | $QSUB -N all_cv_k$K -o all_cv_k$K.stdout -pe mpi 12
done

#when all runs are done get CV scores
grep CV all_cv_maf05_k*.log \
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
    --vcf  ../build_snp_panel/auto_maf_ld.vcf \
    --remove-indv ERR103051 \
    --remove-indv ERR119612 \
    --remove-indv ERR119613 \
    --remove-indv ERR310940 \
    --remove-indv ERR539855 \
    --remove-indv ERR539857 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >haem_auto_maf_ld.vcf

grep "#" haem_auto_maf_ld.vcf \
    | tail -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >haem_auto_maf_ld.list

#convert to ped
plink \
    --vcf haem_auto_maf_ld.vcf \
    --out haem_auto_maf_ld \
    --recode12 \
    --allow-extra-chr

#submit jobs to the cluster
for K in $(seq 1 20); do
    CMD="$WORK_DIR/scripts/admixture_linux-1.3.0/admixture \
        --cv=1000 \
        -j12 \
        haem_auto_maf_ld.ped \
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
# SUPERVISED?
$WORK_DIR/scripts/admixture_linux-1.3.0/admixture \
    -j12 \
    --cv=1000 \
    --supervised \
    haem_auto_maf_ld.ped \
    2 \
    >supervised.log

#plot in R on local computer using the script <R SCRIPT>
