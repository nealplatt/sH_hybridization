#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# build_cohort_vcf - after genotyping, build and filter cohort vcf for 
#       recalibration

# TODO(nplatt): update comments

source master/nplatt/sH_hybridizationscripts/set_env.sh

#rule recal_0_merge_cohort_gvcf:
#    input:
#        {contig}.vcf (expand to all)
#        REFERENCE="data/genome/schHae_v1.fa",
#        rules.get_and_index_genome.output
#    output:
#        protected(COHORT_GVCF="results/recal_0/sH_recal_0_cohort.g.vcf"),
#        temp(GVCF_LIST="results/recal_0/gvcfs.list"),
#        temp(TMP="results/recal_0/tmp_cohort.g.vcf")
#    singularity:
#        "snpCalling_v0.0.5.img"
#    threads:
#        12
#    log:
#        "logs/{rulename}.log"
#    shell:
#        """        
#        ls {input.VCFS}>{output.GVCF_LIST}
#
#        gatk MergeVcfs \
#            -I {output.GVCF_LIST} \
#            -O {output.TMP} \
#            -R {input.REFERENCE}
#
#        gatk SortVcf \
#            -I {output.TMP} \
#            -O {output.COHORT_GVCF}
#        """

#probably have to do serial merging
#
#after merging sort
#
#after sorting - select SNPS
#using this as a guide: https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set


cd /master/nplatt/sH_hybridization
source scripts/set_env.sh

cd $BSRCL_DIR

#-------------------------------------------------------------------------------
# merge all 50 interval vcfs per individual into on vcf per individual

for SAMPLE in $(cat $SAMPLE_LIST); do


SINGULARITY_CMD="/opt/bio/singularity/bin/singularity exec snpCalling_v0.0.4.simg"
QSUB="qsub -V -cwd -j y -S /bin/bash -q all.q"

split -d -n l/1000 --additional-suffix _round0.list gvcfs.list gvcfs.


for i in $(ls gvcfs.*_round0.list); do
  GATK_CMD="gatk --java-options "-Xmx4G" \
    MergeVcfs \
    -I  $i \
    -O round0/$i.g.vcf \
    -R /master/nplatt/sH_hybridization/data/genome/schHae_v1.fa"

  echo $SINGULARITY_CMD $GATK_CMD >$i.sh
        
  cat $i.sh | $QSUB -pe mpi 12 -N $i -o $i.o 
   
done

for LOG in *.o; do
  if [[ $(grep "picard.vcf.MergeVcfs done" $LOG) ]]; then
    PASSED=$((PASSED+1))
  else
    FAILED=$((FAILED+1))
    echo $LOG
    fi 
    TOTAL=$((TOTAL+1))

done

ls round0/*.g.vcf >round2.list

split -d -n l/100 --additional-suffix _round0.list .list round2_gvcfs.

mkdir round1

for i in $(ls  round2_gvcfs.*_round0.list); do
    GATK_CMD="gatk --java-options "-Xmx4G" \
        MergeVcfs \
        -I  $i \
        -O round1/$i.g.vcf \
        -R /master/nplatt/sH_hybridization/data/genome/schHae_v1.fa"

    echo $SINGULARITY_CMD $GATK_CMD >$i.sh
        
    cat $i.sh | $QSUB -pe mpi 12 -N $i -o $i.o 

done

for LOG in *.o; do
    if [[ $(grep "picard.vcf.MergeVcfs done" $LOG) ]]; then
        PASSED=$((PASSED+1))
    else
        FAILED=$((FAILED+1))
        echo $LOG
    fi 
    TOTAL=$((TOTAL+1))

done

rm *list.o *list.sh


gatk --java-options "-Xmx4G" \
        MergeVcfs \
        -I  round2.list \
        -O cohort.g.vcf \
        -R /master/nplatt/sH_hybridization/data/genome/schHae_v1.fa

mv cohort.g.vcf (up)

#now sort
singularity exec snpCalling_v0.0.4.simg gatk --java-options "-Xmx4G" SortVcf -I cohort.g.vcf -O cohort_sorted.g.vcf -TMP_DIR ./tmp_sort

#select all of the SNPs
singularity exec ~/snpCalling_v0.0.5.img gatk SelectVariants -V cohort_sorted.g.vcf -select-type SNP -O cohort_raw_snps.g.vcf -R /master/nplatt/sH_hybridization/data/genome/schHae_v1.fa

#select all of the indels
singularity exec ~/snpCalling_v0.0.5.img gatk SelectVariants -V cohort_sorted.g.vcf -select-type INDEL -O cohort_raw_indels.g.vcf -R /master/nplatt/sH_hybridization/data/genome/schHae_v1.fa

#filter out the low qual SNPs
singularity exec ~/snpCalling_v0.0.5.img gatk VariantFiltration \
    -R -R /master/nplatt/sH_hybridization/data/genome/schHae_v1.fa \
    -V cohort_raw_snps.g.vcf \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "initial_snp_filter" \
    -o filtered_snps.vcf 


