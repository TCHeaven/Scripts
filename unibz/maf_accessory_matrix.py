#!/usr/bin/env python3

import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(
    description="Create genome presence/absence matrix from HAL-derived MAF"
)

parser.add_argument(
    "--maf",
    required=True,
    help="MAF file"
)

parser.add_argument(
    "--window",
    type=int,
    default=5000,
    help="reference window size"
)

parser.add_argument(
    "--reference",
    required=True,
    help="reference genome name in MAF"
)

parser.add_argument(
    "--output",
    default="accessory_matrix.tsv"
)

args = parser.parse_args()


presence = defaultdict(set)

genomes=set()

current_block=[]


def process_block(block):
    if not block:
        return
    ref=None
    for line in block:
        fields=line.split()
        genome=fields[1].split(".")[0]
        genomes.add(genome)
        if genome==args.reference:
            ref=fields
    if ref is None:
        return
    ref_start=int(ref[2])
    ref_len=int(ref[3])
    start_window=ref_start//args.window
    end_window=(ref_start+ref_len)//args.window
    block_genomes=set()
    for line in block:
        fields=line.split()
        genome=fields[1].split(".")[0]
        block_genomes.add(genome)
    for w in range(start_window,end_window+1):
        presence[w].update(block_genomes)

with open(args.maf) as f:
    block=[]
    for line in f:
        if line.startswith("a"):
            process_block(block)
            block=[]
        elif line.startswith("s"):
            block.append(line.strip())
    process_block(block)



genomes=sorted(genomes)


with open(args.output,"w") as out:
    out.write(
        "window\t"
        +"\t".join(genomes)
        +"\n"
    )
    for w in sorted(presence):
        row=[str(w)]
        for g in genomes:
            if g in presence[w]:
                row.append("1")
            else:
                row.append("0")
        out.write(
            "\t".join(row)
            +"\n"
        )


print("Done")
print("Genomes:")
for g in genomes:
    print(g)
