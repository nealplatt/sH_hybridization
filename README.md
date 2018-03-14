# SNP calling in SH exomes.

<explanation/abstract of project>

---
### To Run:

1) Clone repository
```
git clone https://github.com/nealplatt/sH_hybridization.git
```

2) Build singularity image/container
```    
sudo singularity build scripts/Singularity <image name>
```

4) Modify major env variables in the `scripts/set-env.sh` 

3) Download raw sequence reads to data/
```    
# SRA links provided upon publication
```

5) Execute scripts in scripts in the following order:
    ```
    set-env.sh
    map_filterReads.sh
    map_prepGenome.sh
    map_mapAndProcessBam.sh
    map_checkBams.sh
    map_buildListOfIntervals.sh
    bqsr_haplotypeCaller.sh
    bqsr_prepForGenotype.sh
    bqsr_genotype.sh
    bqsr_buildCohortVcf.sh
    bqsr_recalibration.sh
    ```


---

### Dir tree inc. major files
  ```
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
  ```
