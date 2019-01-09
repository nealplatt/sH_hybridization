# clean up headers and trim peptidase fasta file from
# pfam

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

fasta_in = open(sys.argv[1], "r")
fasta_out = open(sys.argv[2], "w")

trimmed_sequences = []  # Setup an empty list

for record in SeqIO.parse(fasta_in, "fasta"):
    peptidase_pos = record.id.split("/")[1]
    peptidase_start = int(peptidase_pos.split("-")[0]) - 1
    peptidase_end = int(peptidase_pos.split("-")[1])
    peptidase_seq = record.seq[peptidase_start:peptidase_end]

    peptidase_seq = SeqRecord(
        peptidase_seq, id=record.id, name=record.name, description=record.description
    )

    trimmed_sequences.append(peptidase_seq)

SeqIO.write(trimmed_sequences, fasta_out, "fasta")
