#...............................................................................
#...............................................................................
#...............................................................................
#ADMIXTURE - needs LD pruned snps and to be in the plink ped format

#NEED TO RE-RUN AND EXCLUDE SEX SNPs

#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

mkdir $RESULTS_DIR/05-admixture
cd $RESULTS_DIR/05-admixture

#convert to ped
plink \
    --vcf $RESULTS_DIR/04-filter/cohort_snps_schMan_final_autosomal_LD.vcf \
    --out cohort_snps_schMan_final_autosomal_LD \
    --recode12 \
    --allow-extra-chr

#UNSUPERVISED - RAW
#The submit 10 jobs to the cluster searching for admixture with K pops
for K in $(seq 1 20); do
    CMD="$WORK_DIR/scripts/admixture_linux-1.3.0/admixture \
        --cv=10 \
        cohort_snps_schMan_final_autosomal_LD.ped \
        $K \
        >cv_auto_k$K.log"

    echo $CMD | $QSUB -N cv_auto_k$K -o cv_auto_k$K.stdout -pe mpi 12
done

#get CV scores
grep CV cv_auto_k*.log \
    | sort -n \
    | awk '{print $3"\t"$4}' \
    | sed 's/(//' \
    | sed 's/)//' \
    | sed 's/://' \
    >admixture_cv_table.tsv
#...............................................................................
#exclude zanzibar samples

vcftools \
    --vcf  $GENO_ALL_DIR/all_schisto_smancoord_filtered_LD-25-5-2.vcf \
    --remove $GENO_ALL_DIR/zanzibar.list \
    --recode \
    --recode-INFO-all \
    --stdout >noTZ_schisto_smancoord_filtered_LD-25-5-2.vcf 

plink \
    --vcf noTZ_schisto_smancoord_filtered_LD-25-5-2.vcf \
    --out noTZ_schisto_smancoord_filtered_LD-25-5-2 \
    --recode12 \
    --allow-extra-chr

#UNSUPERVISED
#The submit 10 jobs to the cluster searching for admixture with K pops
for K in $(seq 1 10); do
    CMD="$WORK_DIR/scripts/admixture_linux-1.3.0/admixture \
        --cv=1000 \
        noTZ_schisto_smancoord_filtered_LD-25-5-2.ped \
        $K \
        >noTZ_cv_k$K.log"

    echo $CMD | $QSUB -N noTZ_cv_k$K -o noTZ_cv_k$K.stdout -pe mpi 12
done



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
