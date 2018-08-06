source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir rehh

cd rehh

#create a list of NE and TZ samples
cp ../niger.list .
cp ../tz.list .
cat niger.list tz.list >samples.list


#polarize vcf by extracting margrebowie vcf
vcftools \
    --vcf ../build_snp_panel/auto_maf.vcf \
    --indv ERR310940 \
    --recode \
    --stdout >marg.vcf

#create a file with the marg sites as 1/1 since these are the only ones
# that need to be changed
grep -v "#" marg.vcf | grep "1/1" | cut -f3 >marg.derived

vcftools \
    --vcf ../beagle/auto_beagle_maf05.vcf \
    --snps marg.derived \
    --recode \
    --stdout \
    >auto_beagle_maf05_marg1.vcf

#use the python script to modify the ancestral alleles.
python ../scripts/python/modify_ancestral_alleles auto_beagle_maf05_marg1.vcf auto_beagle_maf05_marg1_mod0.vcf

#now combine the modified anc alleles with the ones that were already correct
grep -v "#" marg.vcf | grep "0/0" | cut -f3 >marg.anc

vcftools \
    --vcf ../beagle/auto_beagle_maf05.vcf \
    --snps marg.anc \
    --recode \
    --stdout \
    >auto_beagle_maf05_marg0.vcf

#combine the two vcf files and sort
vcfcombine \
    auto_beagle_maf05_marg1_mod0.vcf \
    auto_beagle_maf05_marg0.vcf \
    >anc0.vcf

${ENVIRONMENTS["TITAN SINGULARITY"]} \
     gatk IndexFeatureFile \
        -F anc0.vcf

#add contigs for sman to header
${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk SelectVariants \
        -R $MAN_GENOME \
        -V anc0.vcf \
        -O anc0.vcf_tmp

${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk SortVcf \
        -R $MAN_GENOME \
        -I anc0.vcf_tmp \
        -O auto_beagle_maf05_marg_polarized.vcf
#33,938 sites

#now format the data fro rehh
#change chr names
sed -s 's/SM_V7_//gi' auto_beagle_maf05_marg_polarized.vcf >auto_beagle_maf05_marg_polarized_chr.vcf

for CHR in $(seq 1 7); do
    for POP in  tz niger; do
        
        vcftools \
            --vcf auto_beagle_maf05_marg_polarized_chr.vcf \
            --keep $POP.list \
            --chr $CHR \
            --recode \
            --stdout >$POP"_SM_V7_"$CHR.vcf

        plink \
            --vcf $POP"_SM_V7_"$CHR.vcf \
            --out $POP \
            --recode fastphase 

        grep -v "#" $POP.chr-$CHR.recode.phase.inp \
            | sed -e 's/\(.\)/\1 /g' \
            | sed 1,3d \
            | awk '{print NR" "$0}' \
            >$POP"_"$CHR.inp
    
        grep -v "#" $POP"_SM_V7_"$CHR.vcf \
            | awk '{print $3"\t"$1"\t"$2"\t"$4"\t"$5}' \
            >$POP"_"$CHR.map
    done
done 

#combine into files suitable for whole genome analysis

for CHR in $(seq 1 7); do
    for POP in  tz niger; do
        
       grep -v "#" $POP.chr-$CHR.recode.phase.inp \
            | sed -e 's/\(.\)/\1 /g' \
            | sed 1,3d \
            >$POP"_"$CHR.tmp
    
    done
done 
 
paste niger_?.tmp | awk '{print NR" "$0}' >niger.inp
paste tz_?.tmp | awk '{print NR" "$0}' >tz.inp

cat niger_?.map >niger.map
cat tz_?.map >tz.map
#now analyze in rehhSs

