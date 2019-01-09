import sys

# open necessary files (read in from the cmd line)
# vcf_in=open(sys.argv[1], 'r')
# vcf_out=open(sys.argv[2], 'w')


vcf_in = open("auto_beagle_maf05_marg1.vcf", "r")
vcf_out = open("auto_beagle_maf05_marg1_mod0.vcf", "w")

# ADD THE ANCESTRAL ALLELE TO THE INFO FIELD
for vcf_entry in vcf_in:

    if vcf_entry.startswith("#"):
        # its a header line
        vcf_out.write(vcf_entry)
    else:
        # its an entry and the alt/ref needs to be swapped and 0/1s tr'd
        (chrom, pos, id, ref, alt, qual, filter, info, format) = vcf_entry.split("\t")[
            0:9
        ]

        # get the phased data and replace 0 w/ 1 and 1 w/ 0.
        phased_hap = vcf_entry.split("\t")[9:]
        phased_hap = "\t".join(phased_hap)
        phased_hap = phased_hap.replace("0", "a").replace("1", "0").replace("a", "1")

        # recreate the vcf entry by swaping alt and ref and included the trd phased haps
        vcf_entry = "\t".join(
            [chrom, pos, id, alt, ref, qual, filter, info, format, phased_hap]
        )

        # print
        vcf_out.write(vcf_entry)

# clean up
vcf_in.close()
vcf_out.close()
