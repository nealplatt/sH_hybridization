# VCF TO HAPLOID FASTA
# python vcf_to_haploid_fasta.py <input vcf> <output fasta>
import sys

haplotype_A = {}
haplotype_B = {}


# open necessary files (read in from the cmd line)
vcf_file_in = open(sys.argv[1], "r")

for vcf_entry in vcf_file_in:

    vcf_entry = vcf_entry.rstrip()

    if vcf_entry.startswith("#CHROM"):
        # its the header line

        # create list of samples
        samples = vcf_entry.split("\t")[9:]

        # create dictionary of seqeucnes wtih keys from samples (empty)
        # for haplotype in sequences.fromkeys(samples, ''):
        haplotype_A = haplotype_A.fromkeys(samples, "")
        haplotype_B = haplotype_B.fromkeys(samples, "")

    elif not vcf_entry.startswith("#"):
        ref_allele = vcf_entry.split("\t")[3]
        alt_allele = vcf_entry.split("\t")[4]
        genotypes = vcf_entry.split("\t")[9:]

        for i in range(0, len(genotypes)):
            # print("for loop")
            sample = samples[i]

            # convert haplotype to nucleotide
            if genotypes[i].split("|")[0] == "0":
                hap_A = ref_allele
            else:
                hap_A = alt_allele

            if genotypes[i].split("|")[1] == "0":
                hap_B = ref_allele
            else:
                hap_B = alt_allele

            haplotype_A[sample] += hap_A
            haplotype_B[sample] += hap_B


vcf_file_in.close()

# now pring everything in the sequences dictionary
fasta_file_out = open(sys.argv[2], "w")
for sample in samples:
    fasta_file_out.write(">" + sample + "_A\n" + haplotype_A[sample] + "\n")
    fasta_file_out.write(">" + sample + "_B\n" + haplotype_B[sample] + "\n")


fasta_file_out.close()
