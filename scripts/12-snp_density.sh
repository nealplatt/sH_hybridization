source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

mkdir $RESULTS_DIR/admixture
cd $RESULTS_DIR/admixture

#make windows
bedtools makewindows -g ../../data/genome/Smansoni_v7.fa.fai -w 10000 -s 10000 | grep '\<SM_V7_.\>' >sman_10k-10k.windows

#get vcf files
cp ../11-pop-assign/samples.list .
grep Sh.NE samples.list >niger.samples
grep Sh.TZ samples.list >tz.samples

#all
cp ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf all.vcf

#group
vcftools \
    --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
    --keep samples.list \
    --recode \
    --stdout \
    >haem_group.vcf

#niger
vcftools \
    --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
    --keep samples.list \
    --recode \
    --stdout \
    >niger.vcf

#bovis
vcftools \
    --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
    --indv ERR103048 \
    --recode \
    --stdout \
    >bovis.vcf

#tz
vcftools \
    --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
    --keep samples.list \
    --recode \
    --stdout \
    >tz.vcf

for TYPE in tz bovis niger haem_group all; do
    bcftools view \
        -m2 \
        -M2 \
        -v snps \
        $TYPE.vcf \
        >$TYPE"_biallelic.vcf"

    vcf2bed \
        --do-not-sort \
        --max-mem=2G \
        < $TYPE"_biallelic.vcf" \
        | awk '{print $1"\t"$2"\t"$3"\t"$4}' \
        >$TYPE"_biallelic.bed"
done

#intersect with windows (to get snp/kb counts)
bedtools intersect -c -a sman_10k-10k.windows -b all_biallelic.bed  >sman_10k-10k.counts


#combine into a single file (with a header)


