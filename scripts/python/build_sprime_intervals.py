# build_sprime_intervals.py <sprime score file>
import sys

sprime_score_table_in = open(sys.argv[1], "r")

segment_start = {}
segment_end = {}
segment_chr = {}
snp_ids = {}
score = {}

for introgressed_snp in sprime_score_table_in:

    introgressed_snp = introgressed_snp.rstrip()

    if not introgressed_snp.startswith("CHROM"):
        sman_chr = introgressed_snp.split("\t")[0]
        sman_bp = introgressed_snp.split("\t")[1]
        shae_snp_id = introgressed_snp.split("\t")[2]
        segment = introgressed_snp.split("\t")[5]
        segment_score = introgressed_snp.split("\t")[7]

        if segment not in segment_start:
            segment_start[segment] = sman_bp
            segment_end[segment] = sman_bp
            snp_ids[segment] = [shae_snp_id]
            segment_chr[segment] = sman_chr
            score[segment] = segment_score
        else:
            if segment_start[segment] > sman_bp:
                segment_start[segment] = sman_bp
            if segment_end[segment] < sman_bp:
                segment_end[segment] = sman_bp
            snp_ids[segment].append(shae_snp_id)

sprime_score_table_in.close()

for segment in segment_start:

    # chr start stop name (list of snp ids)
    chrom = segment_chr[segment]
    start = segment_start[segment]
    stop = segment_end[segment]
    segment_id = segment
    snps = snp_ids[segment]
    segment_score = score[segment]
    snp_list = ",".join(snps)

    # print chrom, "\t", start, "\t", stop, "\t", score, "\t", containing_snps
    print(
        chrom,
        "\t",
        start,
        "\t",
        stop,
        "\t",
        segment_id,
        "\t",
        segment_score,
        "\t",
        snp_list,
    )
