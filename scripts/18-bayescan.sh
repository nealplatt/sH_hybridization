source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir bayescan

cd bayescan

cp ../niger.list .
cp ../tz.list .

vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --keep ../niger.list \
    --keep ../tz.list \
    --stdout \
    --recode \
    >tz-ne.vcf
#After filtering, kept 93 out of 96 Individuals
#After filtering, kept 370770 out of a possible 370770 Sites

vcftools \
    --vcf tz-ne.vcf \
    --maf 0.05 \
    --stdout \
    --recode \
    --min-alleles 2 \
    --max-alleles 2 \
    >tz-ne_maf05_bi.vcf
#After filtering, kept 93 out of 93 Individuals
#After filtering, kept 35241 out of a possible 370770 Sites



#converted to bayescan format wtih pgd spider on local computer
#vcf ->pgd ->bayescan  (tz-ne_maf05_bi.bayescan)

wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
unzip BayeScan2.1.zip
chmod u+x BayeScan2.1/binaries/BayeScan2.1_linux64bits 


DATE=$(date +"%Y-%m-%d")

for CHAIN in $(seq -w 1 100); do

    mkdir -p $DATE/chain$CHAIN
    cd $DATE/chain$CHAIN
    
    JID=chain$CHAIN
    LOG=$JID.log
    
    CMD="../../BayeScan2.1/binaries/BayeScan2.1_linux64bits \
        ../../tz-ne_maf05_bi.bayescan \
        -threads 12 \
        -o tz-ne_maf05_bi_bayescan_pr10_pi10k_bin50K_ngen50_npb50_thin20_chain$CHAIN \
        -pr_odds 10 \
        -pilot 10000 \
        -burn 50000 \
        -nbp 50 \
        -n 50000 \
        -thin 20"

    JOB_QSUB=$QSUB" -N $JID -o $LOG -pe mpi 12"

    echo $CMD | $JOB_QSUB 

    cd ../..
done
