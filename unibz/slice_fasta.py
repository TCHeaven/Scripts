#!/usr/bin/env python3
"""
Slice FASTA sequences at regular intervals.

Usage:
    python slice_fasta.py -i genome.fasta -o slices.fasta -s 1000000 -t 6000000 -f 0
"""

import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Slice FASTA sequences at regular intervals")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    parser.add_argument("-s", "--slice-size", type=int, default=1_000_000,
                        help="Length of each slice (default: 1,000,000)")
    parser.add_argument("-t", "--step", type=int, default=6_000_000,
                        help="Distance to next slice start (default: 6,000,000)")
    parser.add_argument("-f", "--offset", type=int, default=0,
                        help="Starting offset in sequence (default: 0)")
    args = parser.parse_args()

    with open(args.output, "w") as out_f:
        for record in SeqIO.parse(args.input, "fasta"):
            seq_len = len(record.seq)
            slice_start = args.offset
            count = 1
            while slice_start < seq_len:
                slice_end = min(slice_start + args.slice_size, seq_len)
                sliced_seq = record.seq[slice_start:slice_end]
                if len(sliced_seq) == args.slice_size:
                    out_f.write(f">{record.id}_slice{count}_{slice_start+1}_{slice_end}\n")
                    out_f.write(str(sliced_seq) + "\n")
                    count += 1
                slice_start += args.step

if __name__ == "__main__":
    main()
