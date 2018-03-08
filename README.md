RNPlatt
6 Feb 2018
calling SNPs in sH with exome capture array data


TODO
1 - bind local directories in singularity image



To Run:
git clone https://github.com/nealplatt/sH_hybridization.git
#download data to data/exome_seq_reads and data/genomes
cd scripts
snakemake --cluster "qsub -pe mpi {threads}" --jobs 200 --latency 240
