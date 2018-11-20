source activate snp_calling
source /master/nplatt/schisto_hybridization/scripts/set_env.sh

cd $RESULTS_DIR/whole_genome_invadolysin_haplotypes

cp ../dnds/Smp_127030_cds_schHae.bed .
cp ../dnds/KL250964.1.fasta .
cp ../dnds/KL250964.1.fasta.fai .
cp ../dnds/KL250964.1.dict .


grep "#" ../build_whole_genome_snp_panel/wg_auto_beagle.vcf \
    | tail -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >whole_genome_samples.list


python ../../scripts/schMan_to_schHae_vcf_coords.py \
    ../build_whole_genome_snp_panel/wg_auto_beagle.vcf \
    wg_auto_beagle_schHae.vcf

python ../../scripts/diploid_to_haploid_vcf.py \
    wg_auto_beagle_schHae.vcf \
    wg_auto_beagle_schHae_phased_hap_A.vcf \
    wg_auto_beagle_schHae_phased_hap_B.vcf


rm whole_genome_Smp_127030_haplotypes_cds.fasta
for SAMPLE in $(cat whole_genome_samples.list); do
    for HAPLOTYPE in A B; do

        #get the vcf
        vcftools \
            --vcf wg_auto_beagle_schHae_phased_hap_$HAPLOTYPE.vcf \
            --indv $SAMPLE \
            --bed ../dnds/Smp_127030_cds_schHae.bed \
            --recode \
            --recode-INFO-all \
            --stdout \
            >$SAMPLE.vcf

        #get the fasta
        java -jar ../../scripts/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar  \
           -T FastaAlternateReferenceMaker \
           -R KL250964.1.fasta \
           -IUPAC $SAMPLE \
           -o $SAMPLE.fasta \
           -V $SAMPLE.vcf

        #clean up fasta for extraction
        sed -i "1c>KL250964.1" $SAMPLE.fasta
        samtools faidx $SAMPLE.fasta

        #extract each cds and them combine into one
        bedtools getfasta \
            -bed Smp_127030_cds_schHae.bed \
            -fi $SAMPLE.fasta \
            -fo $SAMPLE"_exons.fasta"

        #combine into a single sequence
         grep -v "^>" $SAMPLE"_exons.fasta" \
            | awk 'BEGIN { ORS=""; print ">Sequence_name\n" } { print }' > $SAMPLE"_cds.fasta"

        #change the header
         sed -i "1c>$HAPLOTYPE.$SAMPLE" $SAMPLE"_cds.fasta"
        
        #add it to combined fasta file
         cat $SAMPLE"_cds.fasta" >>whole_genome_Smp_127030_haplotypes_cds.fasta
         echo >>whole_genome_Smp_127030_haplotypes_cds.fasta

        #clean up
        rm $SAMPLE.fasta
        rm $SAMPLE.vcf*
        rm $SAMPLE"_exons.fasta"
        rm $SAMPLE"_cds.fasta"
        rm out.log
    done
done
