#!/bin/bash
#
# SNP calling in S. haemotobium hybridzone(s).
# NPlatt
# Neal.platt@gmail.com

# bsrcl_genotype - runs GATK's GenotypeGVCF on each sample/interval combination

# TODO(nplatt): add code to check that all combos run successfully
# TODO(nplatt): update comments


source master/nplatt/sH_hybridization/scripts/set-env.sh

cd $BSRCL_DIR

NUM_JOBS_IN_QUEUE=0
MAX_JOBS_ALLOWED=3000

mkdir interval_vcf

for INTERVAL in $(ls db); do
    GENOTYPE_JOB_NAME=$INTERVAL".genotype"
    THREADS=1

    IN_DB="gendb://db/$INTERVAL"
    OUT_VCF="interval_vcf/$INTERVAL.vcf"
    
    GENOTYPE="$SINGULARITY gatk GenotypeGVCFs -R $REFERENCE -V $IN_DB -new-qual -O $OUT_VCF"

    GENOTYPE_QSUB="$QSUB -pe mpi $THREADS -N $GENOTYPE_JOB_NAME -o logs/$GENOTYPE_JOB_NAME.log"
    echo $GENOTYPE >scripts/$GENOTYPE_JOB_NAME.sh


    NUM_JOBS_IN_QUEUE=$(qstat | grep nplatt | wc -l)

    while [ $NUM_JOBS_IN_QUEUE -gt $MAX_JOBS_ALLOWED ]; do
        sleep 1s
        echo -n "."
        NUM_JOBS_IN_QUEUE=$(qstat | grep nplatt | wc -l)
    done

    cat scripts/$GENOTYPE_JOB_NAME.sh | $GENOTYPE_QSUB
    
done


#################################

#CHECK FOR COMPLETION


