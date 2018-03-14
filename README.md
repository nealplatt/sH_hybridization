#SNP calling in SH exomes.

<explanation/abstract of project>

----------
##To Run:

1) Clone repository

    git clone https://github.com/nealplatt/sH_hybridization.git


2) Download raw sequence reads to data/
    
    # SRA links provided upon publication

3) Modify major env variables in the scripts/set-env.sh 


4) Execute scripts in scripts in the following order:

    set-env.sh
    map_filterReads.sh
    map_prepGenome.sh
    map_mapAndProcessBam.sh
    map_checkBams.sh
    map_buildListOfIntervals.sh
    bqsr-r1_haplotypeCaller.sh
    bqsr-r1_prepForGenotype.sh
    bqsr-r1_genotype.sh
    bqsr-r1_buildCohortVcf.sh
    bqsr-r1_bqsr.sh
    bqsr-r2_haplotypeCaller.sh
    bqsr-r2_prepForGenotype.sh
    bqsr-r2_genotype.sh
    bqsr-r2_buildCohortVcf.sh
    bqsr-r2_bqsr.sh
    bqsr-r3_haplotypeCaller.sh
    bqsr-r3_prepForGenotype.sh
    bqsr-r3_genotype.sh
    # ***bqsr-r3_buildCohortVcf.sh
    # ***bqsr-r3_bqsr.sh



----------
##Dir tree inc. major files
    sH_hybridization
    |----data
    |----scripts
    |--------<see above>
    |----data
    |--------raw_seq_data
    |--------genome
    |----results
    |--------filter_reads
    |--------map_reads
    |--------intervals
    |--------base_recalibration
    |------------r0_bqsr_tables
    |------------r1_bqsr_tables
    |------------rN_bqsr_tables
    |------------<final modified bams>
