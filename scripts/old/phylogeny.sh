#...............................................................................
#...............................................................................
#Phylogeny

#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

mkdir $GENO_ALL_DIR/phylogeny
cd $GENO_ALL_DIR/phylogeny

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



