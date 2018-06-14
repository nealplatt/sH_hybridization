#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir pop_assign

cd pop_assign

mkdir vcf ped
#get vcfs in sliding windows
bedtools makewindows -g ../../data/genome/Smansoni_v7.fa.fai -w 150000 -s 50000 | grep '\<SM_V7_.\>' >sman_150k-50k.windows

#find intervals with atleast 100 snps
vcf2bed \
    --do-not-sort \
    --max-mem=2G \
    < ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
    | awk '{print $1"\t"$2"\t"$3"\t"$4}' \
    >cohort_snps_schMan_autosomal_panel.bed
    #339,022 remaining snps

#intersect all snps with windows and count
bedtools intersect -c -a sman_150k-50k.windows -b cohort_snps_schMan_autosomal_panel.bed | awk '{if($4>=100) print $0}' >sman_150k-50k_n100.windows

#only keep shaem and bovis samples


#window/extract vcfs
while read WINDOW; do
    CHR=$(echo $WINDOW | cut -f1 -d" ") 
    START=$(echo $WINDOW | cut -f2 -d" ")
    STOP=$(echo $WINDOW | cut -f3 -d" ")

    VCF="$CHR:"$START"-"$STOP.vcf
    PED="$CHR:"$START"-"$STOP

    #get vcf
    vcftools \
        --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
        --chr $CHR \
        --from-bp $START \
        --to-bp $STOP \
        --keep samples.list \
        --recode \
        --stdout \
        >vcf/$VCF

    #convert to ped
    plink \
        --vcf vcf/$VCF \
        --out ped/$PED \
        --recode12 \
        --allow-extra-chr

done <sman_150k-50k_n100.windows

#build supervised sample assignment list



#assign pops
for K in $(seq 1 20); do
    CMD="$WORK_DIR/scripts/admixture_linux-1.3.0/admixture \
        --cv=1000 \
        --supervised \
        -j12 \
        <PED_IN> \
        $K \
        >all_cv_k$K.log"

    echo $CMD | $QSUB -N all_cv_k$K -o all_cv_k$K.stdout -pe mpi 12
done

mkdir q_files
mv *.Q q_files

#can't access the files because of the ":" in the file name
cd q_files
rename : - *.Q
