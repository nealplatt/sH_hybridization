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

# remove invariant sites
    raxmlHPC \
        -f a \
        -m ASC_GTRCAT \
        --asc-corr=lewis \
        -p 12345 \
        -x 54321 \
        -# 1000 \
        -s auto_maf_ld.fas \
        -n auto_maf_ld >out.log

#this WILL fail and give a list of invariant sites
cat out.log \
    | grep Site \
    | cut -f2 -d" " \
    | sort -n \
    >invariant.sites

#turn the list of invariant sites into a nice bash command to delete them
######## NOT PRETTY #######
tr '\n' ',' <invariant.sites
cut --complement \
    -c 22,38,42,47,48,57,68,74,79,84,121,137,169,178,181,183,185,186,190,202,217,225,231,243,244,250,251,273,275,289,328,354,371,375,410,411,415,422,429,438,459,472,475,476,479,480,486,513,525,546,555,565,579,614,621,649,667,672,708,710,750,779,792,794,812,824,850,854,859,867,888,889,892,917,921,930,933,938,958,963,971,988,998,1011,1017,1030,1059,1069,1071,1072,1075,1076,1077,1078,1079,1081,1082,1083,1139,1153,1158,1159,1162,1176,1217,1232,1246,1305,1323,1347,1365,1368,1371,1376,1377,1401,1406,1415,1434,1451,1466,1477,1479,1480,1491,1496,1505,1607,1622,1628,1672,1679,1694,1696,1698,1716,1717,1726,1742,1763,1803,1804,1805,1809,1813,1816,1820,1824,1825,1836,1879,1885,1918,1920,1937,1949,1981,2008,2011,2026,2032,2061,2110,2135,2140,2154,2155,2159,2163,2240,2246,2256,2269,2284,2288,2311,2334,2335,2340,2349,2367,2370,2372,2460,2470,2486,2491,2540,2544,2545,2546,2559,2573,2604,2635,2648,2674,2675,2676,2677,2682,2690,2695,2713,2718,2725,2739,2756,2762,2786,2789,2813,2850,2917,2930,2933,2940,2953,2965,2980,3032,3034,3035,3039,3046,3052,3057,3076,3078,3079,3098,3104,3112,3129,3131,3137,3138,3167,3174,3186,3195,3202,3203,3230,3235,3257,3296,3310,3335,3378,3380,3382,3416,3463,3488,3561,3577,3581,3582,3584,3588,3627,3682,3690,3691,3693,3702,3705,3713,3720,3752,3800,3821,3835,3839,3871,3875,3882,3890,3910,4047,4060,4086,4087,4134,4157,4172,4174,4201,4250,4267,4273,4276,4298,4334,4337,4355,4400,4443,4445,4446,4447,4462,4466,4481,4493,4497,4512,4514,4518,4546,4551,4552,4560,4576,4591,4613,4625,4641,4669,4722,4741,4763,4771,4814,4840,4845,4848,4851,4852,4918,4921,4922,4926,4927,4955,4961,4988,5012,5031,5070,5093,5095,5096,5100,5196,5215,5228,5230,5249,5252,5310,5320,5352,5362,5392,5402,5456,5459,5467,5468,5475,5489,5490,5497,5503,5533,5536,5537,5555,5569,5583,5606,5611,5615,5623,5631,5634,5639,5640,5680,5682,5695,5715,5755,5771,5775,5779,5786,5792,5794,5795,5799,5826,5844,5845,5873,5896,5900 \
    auto_maf_ld.fas \
    >auto_maf_ld_var.fas

#submit raxml job
JOB_QSUB=$QSUB" -N raxml_1k -o raxml_1k.stdout -pe mpi 12"

echo "source activate snp_calling; 
    raxmlHPC-PTHREADS \
        -f a \
        -T 12 \
        -m ASC_GTRCAT \
        --asc-corr=lewis \
        -p 12345 \
        -x 54321 \
        -# 1000 \
        -s auto_maf_ld_var.fas \
        -n auto_maf_ld_var" | $JOB_QSUB
#-------------------------------------------------------------------------------

#phylogeny of haplotypes
python $WORK_DIR/scripts/vcf_to_haploid_fasta.py \
    ../beagle/auto_beagle_maf05_LD.vcf \
    auto_beagle_maf05_LD.fas

JOB_QSUB=$QSUB" -N auto_beagle_maf05_LD_raxml.log -o auto_beagle_maf05_LD_raxml.stdout -pe mpi 12"

echo "source activate snp_calling; 
    raxmlHPC-PTHREADS \
        -f a \
        -T 12 \
        -m ASC_GTRCAT \
        --asc-corr=lewis \
        -p 12345 \
        -x 54321 \
        -# 1000 \
        -s auto_beagle_maf05_LD.fas \
        -n auto_beagle_maf05_LD" | $JOB_QSUB

#-------------------------------------------------------------------------------


#sliding window phylogeny
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



