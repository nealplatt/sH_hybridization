####################################################################################################

#to calculate the hybrid index we need snps that differentiate bovis from zanzibar
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir hindex

cd hindex

#get a bovis and zanzibar sample list

grep "#" ../build_snp_panel/cohort_snps_schMan_autosomal_panel.vcf \
    | tail -n 1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    | grep "Sh.TZ" >tz.samples

echo "ERR103048" >bovis.samples 

vcftools \
    --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel_LD-25-5-2.vcf \
    --weir-fst-pop tz.samples \
    --weir-fst-pop bovis.samples \
    --stdout \
    | awk '$3==1 {print $1"\t"$2}' \
    >bovis_tz_fixed_snps_LD.list                
#9266 snps fst=1

#now get only those snps that are fixed between bovis/zanzibar
vcftools \
    --vcf ../build_snp_panel/cohort_snps_schMan_autosomal_panel_LD-25-5-2.vcf \
    --positions bovis_tz_fixed_snps_LD.list  \
    --remove-indv ERR103051 \
    --remove-indv ERR119613 \
    --remove-indv ERR310937 \
    --remove-indv ERR310940 \
    --remove-indv ERR539855 \
    --remove-indv ERR539857 \
    --stdout \
    --recode \
    >bovis_zanzibar_fixed_snps_autosomal_LD.vcf      

#convert this list to table for hindex/introgress
$WORK_DIR/scripts/vcf_to_introgress.py bovis_zanzibar_fixed_snps_autosomal_LD.vcf  bovis_zanzibar_fixed_snps_autosomal_LD.geno bovis_zanzibar_fixed_snps_autosomal_LD.loci

#in the matirx add 2 lines to head with location/ID# and in the loc and loci,type
grep "#" bovis_zanzibar_fixed_snps_autosomal_LD.vcf \
    | tail -n 1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >bovis_zanzibar_fixed_snps_autosomal_LD.samples

#for the sake of speed (meeting tomorrow) did this manually

#reducing data set for speed
 grep -v "#" bovis_zanzibar_fixed_snps_autosomal_LD.vcf  | cut -f3 | shuf | shuf | shuf | shuf | head -n 600 >subsample.list

vcftools \
    --vcf bovis_zanzibar_fixed_snps_autosomal_LD.vcf \
    --snps subsample.list  \
    --stdout \
    --recode \
    >bovis_zanzibar_fixed_snps_autosomal_LD_SUBSAMPLE.vcf 


