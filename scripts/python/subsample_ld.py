from random import random

out_file = open("test.ld", "w")

lines = [
    line for line in open("tz_chr1_snps_schMan_autosomal_panel.ld") if random() >= 0.01
]

for item in lines:
    out_file.write("%s\n" % item)
