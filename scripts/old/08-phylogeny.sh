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
    -c 38,45,46,57,72,76,78,81,118,135,175,178,180,182,183,188,200,214,221,235,236,242,243,276,343,359,392,393,397,404,410,419,424,442,456,457,460,461,469,501,517,537,546,555,568,604,611,641,656,661,693,695,734,762,775,777,794,806,832,837,842,870,871,874,894,898,899,902,911,919,945,953,971,980,993,997,1001,1015,1046,1057,1059,1060,1063,1064,1065,1066,1067,1069,1070,1071,1125,1136,1140,1141,1144,1145,1158,1204,1215,1220,1235,1295,1313,1333,1339,1358,1362,1365,1370,1371,1392,1397,1407,1424,1442,1458,1468,1470,1471,1481,1486,1496,1520,1596,1611,1613,1618,1660,1668,1681,1683,1684,1686,1705,1706,1714,1731,1752,1771,1793,1795,1796,1797,1801,1805,1808,1810,1813,1817,1818,1829,1869,1875,1909,1911,1928,1940,2006,2021,2058,2096,2108,2128,2132,2146,2150,2154,2233,2239,2249,2257,2263,2279,2283,2305,2327,2328,2333,2342,2355,2360,2363,2365,2413,2452,2461,2477,2482,2529,2533,2534,2548,2562,2592,2627,2639,2661,2662,2663,2664,2669,2677,2682,2696,2700,2705,2712,2726,2744,2751,2772,2775,2834,2890,2904,2917,2920,2928,2941,2953,2967,3019,3020,3024,3032,3038,3043,3062,3064,3065,3084,3090,3099,3116,3118,3123,3124,3154,3163,3174,3183,3191,3193,3224,3248,3250,3287,3290,3301,3325,3367,3369,3374,3404,3449,3475,3508,3548,3551,3566,3572,3573,3620,3675,3681,3682,3685,3696,3705,3712,3745,3792,3812,3828,3832,3866,3876,3884,3903,4033,4046,4071,4101,4111,4117,4142,4157,4159,4185,4186,4235,4253,4259,4262,4286,4302,4303,4322,4324,4328,4344,4427,4429,4430,4431,4446,4450,4466,4481,4496,4498,4501,4535,4537,4547,4548,4560,4574,4596,4608,4625,4636,4654,4705,4724,4725,4748,4754,4794,4817,4820,4824,4829,4830,4869,4895,4900,4905,4906,4933,4939,4965,4988,5003,5040,5064,5066,5067,5071,5077,5162,5180,5195,5197,5216,5219,5277,5288,5318,5329,5359,5370,5423,5425,5433,5434,5441,5457,5458,5465,5470,5471,5500,5503,5504,5521,5535,5549,5562,5574,5582,5590,5603,5604,5629,5643,5652,5656,5677,5687,5719,5735,5740,5745,5754,5760,5762,5763,5767,5795,5815,5816,5849,5873,5877 \
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
bedtools makewindows \
    -g ../../data/genome/Smansoni_v7.fa.fai \
    -w 250000 \
    -s 125000 | grep '\<SM_V7_.\>' >sman.windows

#find intervals with atleast 100 snps
vcf2bed \
    --do-not-sort \
    --max-mem=2G \
    < ../beagle/auto_beagle_maf05.vcf \
    | awk '{print $1"\t"$2"\t"$3"\t"$4}' \
    >auto_beagle_maf05.bed
    #339,022 remaining snps

#intersect all snps with windows and count
bedtools intersect \
    -c \
    -a sman.windows \
    -b auto_beagle_maf05.bed \
        | awk '{if($4>=100) print $0}' \
        >sman_n100.windows

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



