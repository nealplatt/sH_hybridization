import sys

# open necessary files (read in from the cmd line)
markers_in = open(sys.argv[1], "r")
vcf_file_in = open(sys.argv[2], "r")
windows_coords_out = open(sys.argv[2], "r")

pos = {}
chr = {}
markers_in = open("control-bov-admixed_maf00.markers.txt", "r")
vcf_file_in = open("../beagle/auto_beagle.vcf", "r")
windows_coords_out = open("test", "w")

SM_V7_1_len = 88881357
SM_V7_2_len = 48130368
SM_V7_3_len = 50458499
SM_V7_4_len = 47279781
SM_V7_5_len = 25256119
SM_V7_6_len = 24989083
SM_V7_7_len = 19288021


chr_length = {}
chr_length["SM_V7_1"] = 0
chr_length["SM_V7_2"] = chr_length["SM_V7_1"] + SM_V7_1_len
chr_length["SM_V7_3"] = chr_length["SM_V7_2"] + SM_V7_2_len
chr_length["SM_V7_4"] = chr_length["SM_V7_3"] + SM_V7_3_len
chr_length["SM_V7_5"] = chr_length["SM_V7_4"] + SM_V7_4_len
chr_length["SM_V7_6"] = chr_length["SM_V7_5"] + SM_V7_5_len
chr_length["SM_V7_7"] = chr_length["SM_V7_6"] + SM_V7_6_len


for vcf_entry in vcf_file_in:

    if not vcf_entry.startswith("#"):
        vcf_info = vcf_entry.split("\t")

        snp_id = vcf_info[2]
        snp_pos = vcf_info[1]
        snp_chr = vcf_info[0]

        pos[snp_id] = snp_pos
        chr[snp_id] = snp_chr

for window_entry in markers_in:
    marker_info = window_entry.split()
    window_name = marker_info[0]
    start_snp = marker_info[1]
    stop_snp = marker_info[-1]
    avg_pos = (int(pos[start_snp]) + int(pos[stop_snp])) / 2
    cul_pos = avg_pos + chr_length[chr[start_snp]]

    entry = "\t".join([window_name, chr[start_snp], str(cul_pos)])

    windows_coords_out.write(entry + "\n")


# clean up
markers_in.close()
vcf_file_in.close()
windows_coords_out.close()
