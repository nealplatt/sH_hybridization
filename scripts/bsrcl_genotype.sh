#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bsrcl_genotype - runs GATK's GenotypeGVCF on each sample/interval combination

# TODO(nplatt): add code to check that all combos run successfully
# TODO(nplatt): update comments

source master/nplatt/sH_hybridizationscripts/set_env.sh

cd $BSRCL_DIR

mkdir interval_vcf

for SAMPLE in $(cat $SAMPLE_LIST); do
    for INTERVAL in $(ls db); do
        GENOTYPE_JOB_NAME=$INTERVAL
        THREADS=1

        IN_DB="gendb://db/$INTERVAL"
        OUT_VCF="interval_vcf/$INTERVAL.vcf"
    
        GENOTYPE="$SINGULARITY gatk GenotypeGVCFs -R $REFERENCE -V $IN_DB -new-qual -O $OUT_VCF"

        $GENOTYPE_QSUB="$QSUB -pe mpi $THREADS -N $GENOTYPE_JOB_NAME -o logs/$GENOTYPE_JOB_NAME.log"
        echo $GENOTYPE >scripts/$GENOTYPE_JOB_NAME.sh

        cat scripts/$GENOTYPE.sh | $GENOTYPE
    done
done

#################################

#CHECK FOR COMPLETION


