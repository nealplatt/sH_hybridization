import sys

# open necessary files (read in from the cmd line)
vcf_in = open(sys.argv[1], "r")
snps_out = open(sys.argv[2], "w")
# vcf_in=open("auto_maf.vcf", 'r')
# snps_out=open("autap.snps", 'w')

informative_snps = {}

for vcf_entry in vcf_in:

    if not vcf_entry.startswith("#"):
        vcf_entry = vcf_entry.rstrip()
        # its an entry and the alt/ref needs to be swapped and 0/1s tr'd
        (chrom, pos, id, ref, alt, qual, filter, info, format) = vcf_entry.split("\t")[
            0:9
        ]
        gt = vcf_entry.split("\t")[9:]

        # get the alleles for each "pop"
        egy_gts = [i.split(":")[0] for i in gt][0]
        bov_gts = [i.split(":")[0] for i in gt][1]
        gui_gts = [i.split(":")[0] for i in gt][3]
        int_gts = [i.split(":")[0] for i in gt][4]
        cur_gts = [i.split(":")[0] for i in gt][5]
        mar_gts = [i.split(":")[0] for i in gt][6]
        ne_gts = [i.split(":")[0] for i in gt][9:54]
        tz_gts = [i.split(":")[0] for i in gt][55:101]
        mat_gts = [
            [i.split(":")[0] for i in gt][2],
            [i.split(":")[0] for i in gt][7],
            [i.split(":")[0] for i in gt][8],
        ]

        # make a composite outgroup set of gts
        out_gts = [gui_gts, int_gts, mar_gts] + mat_gts

        # identify alleles for each locus
        bov_alleles = list(set(bov_gts.split("/")))
        cur_alleles = list(set(cur_gts.split("/")))
        tz_alleles = list(set([y for x in tz_gts for y in x.split("/")]))
        ne_alleles = list(set([y for x in ne_gts for y in x.split("/")]))
        out_alleles = list(set([y for x in out_gts for y in x.split("/")]))

        # now find out if allele is bovis specific (or cur specific)
        for allele in bov_alleles:
            if (
                allele not in cur_alleles
                or allele not in out_alleles
                and allele is not "."
            ):
                snps_out.write(id + "\tbov\t" + allele + "\n")
                # snps_out.write(id + "bov\t" + allele + "\t" + cur_alleles+ "\t" + out_alleles +'\n')
        for allele in cur_alleles:
            if (
                allele not in bov_alleles
                or allele not in out_alleles
                and allele is not "."
            ):
                snps_out.write(id + "\tcur\t" + allele + "\n")
                # snps_out.write(id + "cur\t" + allele + "\t" + bov_alleles+ "\t" + out_alleles +'\n')

vcf_in.close()
snps_out.close()
