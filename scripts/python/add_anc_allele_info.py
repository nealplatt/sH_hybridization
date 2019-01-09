import sys

# open necessary files (read in from the cmd line)
site_in = open(sys.argv[1], "r")
vcf_in = open(sys.argv[2], "r")
vcf_out = open(sys.argv[3], "w")


# vcf_in=open("cohort_snps_schMan_final_autosomal.vcf", 'r')
# site_in=open("anc.alleles", 'r')
# vcf_out=open("aa.vcf", 'w')

anc_allele = {}

for site_entry in site_in:

    chr = site_entry.split("\t")[0].rstrip()
    pos = site_entry.split("\t")[1].rstrip()
    allele = site_entry.split("\t")[2].rstrip()

    name = ":".join([chr, pos])
    anc_allele[name] = allele

site_in.close()

# ADD THE ANCESTRAL ALLELE TO THE INFO FIELD
for vcf_entry in vcf_in:

    if vcf_entry.startswith("#"):
        # its a header line
        # print("header")
        vcf_out.write(vcf_entry)
    else:
        # print("vcf_entry")
        chr = vcf_entry.split("\t")[0]
        pos = vcf_entry.split("\t")[1]
        name = ":".join([chr, pos])

        if name in anc_allele:
            # print('anc_allele')
            aa = "".join(["AA=", anc_allele[name]])
            vcf_front = "\t".join(vcf_entry.split("\t")[:7])
            vcf_back = "\t".join(vcf_entry.split("\t")[7:])

            vcf_back_aa = ";".join([aa, vcf_back])
            vcf_entry = "\t".join([vcf_front, vcf_back_aa])

        # print(vcf_entry_post_mod)
        vcf_out.write(vcf_entry)


# clean up
vcf_in.close()
vcf_out.close()
