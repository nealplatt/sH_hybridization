#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N mpi_raxml
#$ -o mpi_raxml.o$JOB_ID                                                                                                           
#$ -j y
#$ -q all.q
#$ -pe mpi 12


mpirun -np 144 ~/bin/standard-RAxML/raxmlHPC-MPI \
    -m GTRCAT -V \
    -p 12345 \
    -s test_snps.fas \
    -n MPI_test_snps_infer_ML \
    -f a\
    -# 1000 \
    -T 36
