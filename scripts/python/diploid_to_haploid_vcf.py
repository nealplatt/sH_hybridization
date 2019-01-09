import sys

# open necessary files (read in from the cmd line)
vcf_in = open(sys.argv[1], "r")
hap_A_vcf_out = open(sys.argv[2], "w")
hap_B_vcf_out = open(sys.argv[3], "w")
# vcf_in=open("auto_beagle_all_samples.vcf", 'r')
# hap_A_vcf_out=open("hap_A.vcf", 'w')
# hap_B_vcf_out=open("hap_B.vcf", 'w')


for vcf_entry in vcf_in:
    if vcf_entry.startswith("#"):
        hap_A_vcf_out.write(vcf_entry)
        hap_B_vcf_out.write(vcf_entry)

    if not vcf_entry.startswith("#"):
        vcf_entry = vcf_entry.rstrip()

        # split the gts and locus info
        locus_info = vcf_entry.split("\t")[0:9]
        genotypes = vcf_entry.split("\t")[9:]

        # take the gts and split them into an A and B haplotype
        haplotype_A = []
        haplotype_B = []

        for genotype in genotypes:
            (allele_A, allele_B) = genotype.split("|")
            fake_homozygous_A = allele_A + "|" + allele_A
            fake_homozygous_B = allele_B + "|" + allele_B

            haplotype_A.append(fake_homozygous_A)
            haplotype_B.append(fake_homozygous_B)

        # now print each new "haplotype" to the vcf
        new_haplotype_A_entry = "\t".join(locus_info + haplotype_A + ["\n"])
        new_haplotype_B_entry = "\t".join(locus_info + haplotype_B + ["\n"])

        hap_A_vcf_out.write(new_haplotype_A_entry)
        hap_B_vcf_out.write(new_haplotype_B_entry)


vcf_in.close()
hap_A_vcf_out.close()
hap_B_vcf_out.close()
