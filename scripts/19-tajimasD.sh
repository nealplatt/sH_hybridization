#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir diversity

cd diversity

cp ../niger.list .
cp ../tz.list .

vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --stdout \
    --TajimaD 25000 \
    --keep niger.list \
    | awk '{if ($3>=10) print $0}' \
    >ne_w25k_maf00_tajd.txt

vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --stdout \
    --TajimaD 25000 \
    --keep tz.list \
    | awk '{if ($3>=10) print $0}' \
    >tz_w25k_maf00_tajd.txt



vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --stdout \
    --TajimaD 25000 \
    --maf 0.05 \
    --keep niger.list \
    | awk '{if ($3>=10) print $0}' \
    >ne_w25k_maf05_tajd.txt

vcftools \
    --vcf ../beagle/auto_beagle.vcf \
    --stdout \
    --TajimaD 25000 \
    --maf 0.05 \
    --keep tz.list \
    | awk '{if ($3>=10) print $0}' \
    >tz_w25k_maf05_tajd.txt
