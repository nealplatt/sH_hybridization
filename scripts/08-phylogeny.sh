#...............................................................................
#...............................................................................
#Phylogeny

#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

mkdir $RESULTS_DIR/phylogeny
cd $RESULTS_DIR/phylogeny

#convert vcf to fasta
python $WORK_DIR/scripts/vcf_to_diploid_fasta.py \
    ../build_snp_panel/auto_maf_ld.vcf \
    auto_maf_ld.fas

JOB_QSUB=$QSUB" -N raxml_1k -o raxml_1k.stdout -pe mpi 12"

echo "raxmlHPC \
    -f a \
    -T 12 \
    -m GTRCAT \
    -p 12345 \
    -x 54321 \
    -# 1000 \
    -s auto_maf_ld.fas \
    -n auto_maf_ld" | $JOB_QSUB


#-------------------------------------------------------------------------------


JOB_QSUB=$QSUB" -N raxml_1k -o raxml_1k.stdout -pe mpi 12"

echo "raxmlHPC \
    -f a \
    -T 12 \
    -m GTRCAT \
    -p 12345 \
    -x 54321 \
    -# 1000 \
    -s cohort_snps_schMan_autosomal_panel_LD-25-5-2.fas \
    -n cohort_snps_schMan_autosomal_panel_LD-25-5-2_1k" | $JOB_QSUB



#bedtools make windows
bedtools makewindows -g ../../data/genome/Smansoni_v7.fa.fai -w 150000 -s 75000 | grep '\<SM_V7_.\>' >sman.windows

#find intervals with atleast 100 snps
vcf2bed \
    --do-not-sort \
    --max-mem=2G \
    < ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
    | awk '{print $1"\t"$2"\t"$3"\t"$4}' \
    >cohort_snps_schMan_autosomal_panel.bed
    #339,022 remaining snps

#intersect all snps with windows and count
bedtools intersect -c -a sman.windows -b cohort_snps_schMan_autosomal_panel.bed | awk '{if($4>=100) print $0}' >sman_n100.windows

mkdir vcf
#window/extract vcfs
while read WINDOW; do
    CHR=$(echo $WINDOW | cut -f1 -d" ") 
    START=$(echo $WINDOW | cut -f2 -d" ")
    STOP=$(echo $WINDOW | cut -f3 -d" ")

    FILENAME="$CHR:"$START"-"$STOP.vcf

    vcftools \
        --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
        --chr $CHR \
        --from-bp $START\
        --to-bp $STOP\
        --recode \
        --stdout \
        >vcf/$FILENAME
done <sman_n100.windows

#vcf to fasta
mkdir fasta

for VCF_FILE in $(ls vcf/SM_V7*.vcf); do
    
    FASTA_FILE=fasta/$(basename $VCF_FILE .vcf).fas

    python $WORK_DIR/scripts/vcf_to_diploid_fasta.py \
        $VCF_FILE \
       $FASTA_FILE

done 


for FASTA_FILE in $(ls fasta/*.fas); do
    
    RAXML_BASE=$(basename $FASTA_FILE .fas)

    echo "raxmlHPC \
        -f a \
        -m GTRCAT \
        -p 12345 \
        -x 54321 \
        -# 1000 \
        -s $FASTA_FILE \
        -n $RAXML_BASE" | $QSUB

done 

#clean up the directory
mkdir sliding_window
mv RAxML_.SM_V7* sliding_window/





awk '{sum+=$2} END {print sum/NR}' sites_150kb.list 
#166.8 avg sites/window

#raxml tree
#check tree topology






#phylogeny
#doesn't seem like its going to run....trying examl

##convert to phylip - manually and then add smansoni as outgroup
#vcf2bed --do-not-sort <$GENO_ALL_DIR/all_schisto_smancoord_filtered_LD-25-5-2.vcf \
#    | awk '{print $1":"$3"-"$3}' >mansoni_LD_filtered.list

#xargs samtools faidx $MAN_GENOME <mansoni_LD_filtered.list \
#    | grep -v ">" >mansoni_LD_filtered.nucs

#paste -s -d',' mansoni_LD_filtered.nucs \
#    | sed 's/,//g' >mansoni_LD_filtered.fas

#sed -i 's/^/schMan_v7\t/' mansoni_LD_filtered.fas


#cat all_schisto_smancoord_filtered_LD-25-5-2.phy \
#    mansoni_LD_filtered.fas \
#    >all_schisto_smancoord_filtered_LD-25-5-2_smanout.phy

#nano all_schisto_smancoord_filtered_LD-25-5-2_smanout.phy

python $WORK_DIR/scripts/vcf_to_diploid_fasta.py \
    $GENO_ALL_DIR/all_schisto_smancoord_filtered.vcf \
    all_schisto_smancoord_filtered.fas

#convert to phylip (manually)
paste -d "\t" - - <all_schisto_smancoord_filtered.fas \
    >all_schisto_smancoord_filtered.phy

SEQ=$(head -n 1 all_schisto_smancoord_filtered.phy  | cut -f2)                                                                     
echo ${#SEQ} #95855

sed -i '1s/^/ 94 95855\n/' all_schisto_smancoord_filtered.phy

#build binary alignment for examl
$WORK_DIR/scripts/parse-examl \
    -s all_schisto_smancoord_filtered.phy \
    -n all_schisto_smancoord_filtered \
    -m DNA




#build 50 randomized stepwise addition order parsimony trees and 50 random trees
for i in $(seq -w 1 50); do
#    raxmlHPC \
#        -y \
#        -d \
#        -m GTRCAT \
#        -p $RANDOM \
#        -s all_schisto_smancoord_filtered.phy \
#        -n pars_$i

#    raxmlHPC \
#        -y \
#        -d \
#        -m GTRCAT \
#        -p $RANDOM \
#        -s all_schisto_smancoord_filtered.phy \
#        -n rand_$i

    PARS_CMD="mpirun -np 30 /master/nplatt/schisto_hybridization/scripts/examl \
        -t RAxML_randomTree.rand_$i \
        -s all_schisto_smancoord_filtered.binary \
        -n rand_$i \
        -m GAMMA"

    RAND_CMD="mpirun -np 30 /master/nplatt/schisto_hybridization/scripts/examl \
        -t RAxML_randomTree.pars_$i \
        -s all_schisto_smancoord_filtered.binary \
        -n pars_$i \
        -m GAMMA"

    #limit num running jobs
    NUM_JOBS=$(qstat | grep nplatt | wc -l)

    while [ $NUM_JOBS -gt 100 ]; do
        echo "."
        sleep 1
        NUM_JOBS=$(qstat | grep nplatt | wc -l)
    done

    echo PARS_CMD | $QSUB" -N pars_$i -o pars_$i.stdout -pe mpi 30"
    echo RAND_CMD | $QSUB" -N rand_$i -o rand_$i.stdout -pe mpi 30"

done

#from these trees find the best tree
#build examl tree from each --then find best in raxml

#generate 100 alignments for bootstrapping
raxmlHPC \
    -#100 \
    -b 12345 \
    -f j \
    -m GTRCAT \
    -s all_schisto_smancoord_filtered.phy \
    -n REPS

for i in $(seq 0 99); do
    #create the random starting tree (BS)
    raxmlHPC \
        -y \
        -d \
        -m GTRCAT \
        -p $RANDOM \
        -s RAxML_randomTree.bsrep_$i \
        -n bsrep_$i

    #create the examl binary alignment (BS)
    $WORK_DIR/scripts/parse-examl \
        -s all_schisto_smancoord_filtered.phy.BS$i \
        -n all_schisto_smancoord_filtered.phy.BS$i \
        -m DNA

    #Find the tree (BS)
    BS_CMD="mpirun -np 30 /master/nplatt/schisto_hybridization/scripts/examl \
        -t RAxML_randomTree.bsrep_$i \
        -s all_schisto_smancoord_filtered.phy.B"$i".binary \
        -n boot_$i \
        -m GAMMA"
    
    NUM_JOBS=$(qstat | grep nplatt | wc -l)

    while [ $NUM_JOBS -gt 100 ]; do
        echo "."
        sleep 1
        NUM_JOBS=$(qstat | grep nplatt | wc -l)
    done

    echo BS_CMD | $QSUB" -N bs_$i -o bs_$i.stdout -pe mpi 30"
done

#once we have the bs trees we can quantify bipartition frequencies on the best
# tree

vcftools \
    --max-alleles 2 \
    --vcf ../autosomes_schisto_cohort_filtered_LD.vcf  \
    --recode \
    --recode-INFO-all \
    --stdout \
    >autosomes_schisto_cohort_filtered_LD_biallelic.vcf 

python $WORK_DIR/scripts/vcf_to_diploid_fasta.py \
    autosomes_schisto_cohort_filtered_LD_biallelic.vcf \
    autosomes_schisto_cohort_filtered_LD_biallelic.fas

raxmlHPC-PTHREADS \
    -T 10 \
    -f a \
    -m GTRCAT \
    -p 12345 \
    -x 12345 \
    -# 10000 \
    -s autosomes_schisto_cohort_filtered_LD_biallelic.fas \
    -n rapidBoot_fullDataset 
mkdir treeSearch
mkdir bootstrap

cd treeSearch

raxmlHPC \
    -T 12 \
    -m GTRCAT \
    -p 12345 \
    -s ../autosomes_schisto_cohort_filtered_LD_biallelic.fas \
    -# 1000 \
    -n  schistoSNP

cd ../bootstrap

raxmlHPC \
    -T 12 \
    -m GTRCAT \
    -p 12345 \
    -b 12345 \
    -s ../autosomes_schisto_cohort_filtered_LD_biallelic.fas \
    -# 10000 \
    -n  schistoSNP



