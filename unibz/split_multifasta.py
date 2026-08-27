#!/usr/bin/env python3

from pathlib import Path
import re
import sys

input_fasta = Path(sys.argv[1])
output_dir = Path(sys.argv[2])
output_dir.mkdir(parents=True, exist_ok=True)


def clean_filename(header):
    name = header.strip()
    name = re.sub(r'[<>:"/\\|?*\x00-\x1f]', '_', name)
    name = re.sub(r'\s+', '_', name)
    return name[:240] + ".fasta"


header = None
sequence = []

with input_fasta.open() as infile:
    for line in infile:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            # Write previous sequence
            if header is not None:
                outfile = output_dir / clean_filename(header)
                with outfile.open("w") as out:
                    out.write(f">{header}\n")
                    seq = "".join(sequence)
                    for i in range(0, len(seq), 80):
                        out.write(seq[i:i+80] + "\n")
            header = line[1:]
            sequence = []
        else:
            sequence.append(line)


# Write final sequence
if header is not None:
    outfile = output_dir / clean_filename(header)
    with outfile.open("w") as out:
        out.write(f">{header}\n")
        seq = "".join(sequence)
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")


print(f"Done. FASTA files written to: {output_dir}")
