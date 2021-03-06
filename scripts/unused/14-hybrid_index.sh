####################################################################################################

#to calculate the hybrid index we need snps that differentiate bovis from zanzibar
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir hindex

cd hindex

#get a bovis and zanzibar sample list

grep "#" ../build_snp_panel/auto_maf_ld.vcf \
    | tail -n 1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    | grep "Sh.TZ" >tz.samples

echo "ERR103048" >bovis.samples 

vcftools \
    --vcf ../build_snp_panel/auto_maf_ld.vcf  \
    --weir-fst-pop tz.samples \
    --weir-fst-pop bovis.samples \
    --stdout \
    | awk '$3==1 {print $1"\t"$2}' \
    >bovis_tz_fixed_snps_LD.list                
#643 snps fst=1

#now get only those snps that are fixed between bovis/zanzibar
vcftools \
    --vcf ../build_snp_panel/auto_maf_ld.vcf \
    --positions bovis_tz_fixed_snps_LD.list  \
    --remove-indv ERR084970 \
    --remove-indv ERR103051 \
    --remove-indv ERR119612 \
    --remove-indv ERR119613 \
    --remove-indv ERR310937 \
    --remove-indv ERR310940 \
    --remove-indv ERR539855 \
    --remove-indv ERR539857 \
    --remove-indv SRR433865 \
    --stdout \
    --recode \
    >bovis_zanzibar_fixed_snps_autosomal_LD.vcf      

#convert this list to table for hindex/introgress
python $WORK_DIR/scripts/vcf_to_introgress.py bovis_zanzibar_fixed_snps_autosomal_LD.vcf  bovis_zanzibar_fixed_snps_autosomal_LD.geno bovis_zanzibar_fixed_snps_autosomal_LD.loci

#in the matirx add 2 lines to head with location/ID# and in the loc and loci,type
grep "#" bovis_zanzibar_fixed_snps_autosomal_LD.vcf \
    | tail -n 1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >bovis_zanzibar_fixed_snps_autosomal_LD.samples

#for the sake of speed (meeting tomorrow) did this manually



