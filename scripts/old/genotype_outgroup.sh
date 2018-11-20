cd /master/nplatt/schisto_hybridization/results/04-genotype_outgroups/

QSUB="qsub -V -cwd -S /bin/bash -q all.q"    

for SAMPLE in ERR084970 ERR037800 SRR433865 Sm.BR_1278.1 Sm.BR_0447.1 Sm.BR_2039.1 ERR119622 ERR103048 ERR119623; do

    BAM="/master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/"$SAMPLE"_processed.bam"

    N=$(echo $SAMPLE | sed s'/\./_/gi')
    O=$N.stdout
    E=$N.stderr
    THREADS=12
    
    BAM_PATH="/master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/"
  
    SAMPLE_QSUB=$QSUB" -N $N -o $O -e $E -pe mpi $THREADS"
    REFERENCE="/master/nplatt/schisto_hybridization/data/genome/schHae_v1.fa"
    LIST="/master/nplatt/schisto_hybridization/results/04-build_refpanel/sH_ref_panel_snps.list"
    IMAGE="/master/nplatt/schisto_hybridization/config/snpCalling_v0.0.8.img"

    CMD="singularity exec $IMAGE gatk HaplotypeCaller -I $BAM  -O $SAMPLE.vcf -R $REFERENCE -L $LIST"

    echo $CMD | $SAMPLE_QSUB

done


CMD="singularity exec /master/nplatt/schisto_hybridization/config/snpCalling_v0.0.8.img gatk HaplotypeCaller -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//ERR037800_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//ERR084970_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//SRR433865_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//Sm.BR_1278.1_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//Sm.BR_0447.1_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//Sm.BR_2039.1_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//ERR119622_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//ERR103048_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//ERR119623_processed.bam -O outgoups_vs_panel.vcf -R /master/nplatt/schisto_hybridization/data/genome/schHae_v1.fa -L /master/nplatt/schisto_hybridization/results/04-build_refpanel/sH_ref_panel_snps.list"

echo $CMD | $QSUB -N outgoups_vs_panel -o outgoups_vs_panel.stdout -e outgoups_vs_panel.stderr -pe mpi 12

singularity exec /master/nplatt/schisto_hybridization/config/snpCalling_v0.0.8.img \
    gatk CollectReadCounts \
        -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/ERR037800_processed.bam \
        -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/ERR084970_processed.bam \
        -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/SRR433865_processed.bam \
        -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/Sm.BR_1278.1_processed.bam \
        -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/Sm.BR_0447.1_processed.bam \
        -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/Sm.BR_2039.1_processed.bam \
        -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/ERR119622_processed.bam \
        -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/ERR103048_processed.bam \
        -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//ERR119623_processed.bam \
        -O outgoups_depth.vcf \
        --format TSV \
        -L /master/nplatt/schisto_hybridization/results/04-build_refpanel/sH_ref_panel_snps.list

singularity exec /master/nplatt/schisto_hybridization/config/snpCalling_v0.0.8.img gatk HaplotypeCaller -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//ERR037800_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//ERR084970_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//SRR433865_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//Sm.BR_1278.1_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//Sm.BR_0447.1_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//Sm.BR_2039.1_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//ERR119622_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//ERR103048_processed.bam -I /master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//ERR119623_processed.bam -O outgoups_vs_panel_min4.vcf -R /master/nplatt/schisto_hybridization/data/genome/schHae_v1.fa -L /master/nplatt/schisto_hybridization/results/04-build_refpanel/sH_ref_panel_snps.list --minReadsPerAlignmentStart 4 


HAE_SRA_SAMPLES = ["ERR084970", "ERR037800", "SRR433865"]
MAN_SAMPLES     = ["Sm.BR_1278.1", "Sm.BR_0447.1", "Sm.BR_2039.1"]
BOV_SAMPLES     = ["ERR119622", "ERR103048"]
CUR_SAMPLES     = ["ERR119623"]



        singularity exec snpCalling_v0.0.8.img \
            gatk --java-options \"-Xmx4g -Xms4g\" GenomicsDBImport \
                -V {input.VCF_LIST} \
                --genomicsdb-workspace-path {output} \
                -L {wildcards.contig} \
                --reader-threads {threads} \
                --batch-size 24



/master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/ERR037800_processed.bam \
/master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/ERR084970_processed.bam \
/master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/SRR433865_processed.bam \
/master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/Sm.BR_1278.1_processed.bam \
/master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/Sm.BR_0447.1_processed.bam \
/master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/Sm.BR_2039.1_processed.bam \
/master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/ERR119622_processed.bam \
/master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads/ERR103048_processed.bam \
/master/nplatt/schisto_hybridization/results/01-processed_reads/sh_mapped_reads//ERR119623_processed.bam \

