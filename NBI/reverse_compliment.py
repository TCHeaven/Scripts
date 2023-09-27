import sys
import os
import re
from Bio import SeqIO

input_file = sys.argv[1]
output_file = sys.argv[2]

def reverse_complement_fasta(input_file, output_file):
    with open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            record.seq = record.seq.reverse_complement()
            SeqIO.write(record, output_handle, "fasta")

reverse_complement_fasta(input_file, output_file)
