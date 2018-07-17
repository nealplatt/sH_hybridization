#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir ld

cd ld

mkdir schMan_ld
for i in $(seq 1 7); do
    CHR=SM_V7_"$i"

        #get vcf chromosome
        vcftools \
            --vcf ../build_snp_panel/auto_maf.vcf \
            --keep haem.list \
            --maf 0.1 \
            --chr $CHR \
            --recode \
            --stdout \
            >$CHR.vcf

        #convert to ped
        plink \
            --threads 6 \
            --vcf $CHR.vcf \
            --out $CHR \
            --recode12 \
            --allow-extra-chr

        #calc R2 for all snps on the chr
        plink \
            --threads 6 \
            --r2 \
            --file $CHR \
            --out $CHR \
            --allow-extra-chr \
            --ld-window-r2 0 \
            --ld-window 75000 \
            --ld-window-kb 100000


        awk '{print $0"\t"$5-$2}' $CHR.ld >schMan_ld/$CHR.ld
    rm $CHR.ped $CHR.map $CHR.nosex $CHR.log $CHR.vcf $CHR.ld
done

#build new vcf file with schHaem coords
grep -v "#" ../build_snp_panel/auto_maf.vcf \
    | awk  '{print $3"\t"$0}' \
    | sed 's/:/\t/' \
    | awk  '{$3=$3":"$4;; print $0}' \
    | sed 's/ /\t/g' \
    | cut --complement -f4,5 \
    >vcf.entries

grep "#" ../build_snp_panel/auto_maf.vcf | grep -v "contig=<ID=" >vcf.sm_header   

cat vcf.sm_header vcf.entries >broke_header.vcf

#add contigs for sman to header
${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk SelectVariants \
        -R $HAE_GENOME \
        -V broke_header.vcf \
        -O header.vcf
        #75,704 loci remaining

#now sort
${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk SortVcf \
        -R $HAE_GENOME \
        -I header.vcf \
        -O auto_maf_schHae.vcf
        #75,704 snp loci across all (sman homologous) chromosomes

${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk SortVcf \
        -R $HAE_GENOME \
        -I header.vcf \
        -O auto_maf_schHae.vcf
#create list of contigs with snps

bcftools view \
    -m2 \
    -M2 \
    -v snps \
    auto_maf_schHae.vcf \
    >auto_maf_schHae_bi.vcf

vcftools \
    --vcf auto_maf_schHae.vcf \
    --keep haem.list \
    --maf 0.1 \
    --recode \
    --stdout \
    >auto_maf_schHae_bi_haem.vcf

grep -v "#" auto_maf_schHae_bi_haem.vcf \
    | cut -f1 \
    | sort \
    | uniq -c \
    | awk '{ if ($1>1) print $2}' \
    >schHae_contigs.list                                                   

mkdir schHae_ld
#calc ld for all schHae contigs
for CHR in $(cat schHae_contigs.list); do

    #get vcf chromosome
    vcftools \
        --vcf auto_maf_schHae.vcf \
        --keep haem.list \
        --maf 0.1 \
        --chr $CHR \
        --recode \
        --stdout \
        >$CHR".vcf"

    #convert to ped
    plink \
        --threads 6 \
        --vcf $CHR".vcf" \
        --out $CHR \
        --recode12 \
        --allow-extra-chr

    #calc R2 for all snps on the chr
    plink \
        --threads 6 \
        --r2 \
        --file $CHR \
        --out $CHR \
        --allow-extra-chr \
        --ld-window-r2 0 \
        --ld-window 75000 \
        --ld-window-kb 100000

    #calc distance
    awk '{print $0"\t"$5-$2}' $CHR.ld >schHae_ld/$CHR.ld
    #rm $CHR.ped $CHR.map $CHR.nosex $CHR.log $CHR.vcf $CHR.ld

done

rm *temp*

#when finished
cat schHae_ld/*.ld | awk '{print $7","$8}' | grep -v "R2" >schHae_ld.csv
tar -czf schHae_ld.tgz schHae_ld/
rm -r schHae_ld/
#produces a rather small ~75MB file

for i in $(seq 1 7); do

    CHR=SM_V7_"$i"

    awk '{print $7","$8 }' $CHR"_haem.ld" | grep -v R2 >$CHR"_haem.ld.csv" &

done
wait




