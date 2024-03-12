#!/usr/bin/env python3

import sys
import argparse
from collections import defaultdict

ap = argparse.ArgumentParser(description="Extract annotations for genes and build Fischer's contingency tables.")
ap.add_argument('--interpro', required=True, type=str, help='Interproscan annotations in TSV')
ap.add_argument('--set1_name', required=True, type=str, help='Name for genes in set 1')
ap.add_argument('--set1', required=True, type=str, help='File containing list of genes in set 1')
ap.add_argument('--set2_name', required=True, type=str, help='Name for genes in set 2')
ap.add_argument('--set2', required=True, type=str, help='File containing list of genes in set 2')
ap.add_argument('--outdir', required=True, type=str, help='Output directory for results')

conf = ap.parse_args()

Set1_name = conf.set1_name
Set2_name = conf.set2_name

with open(conf.interpro) as f:
    interpro_lines = f.readlines()

with open(conf.set1, 'r') as f:
    set1_genes = set(f.read().splitlines())

with open(conf.set2, 'r') as f:
    set2_genes = set(f.read().splitlines())

print(conf.interpro)
print(Set1_name)
#print(set1_genes)
print(Set2_name)
#print(set2_genes)

set1_count = 0
set2_count = 0
set1_dict = defaultdict(int)
set2_dict = defaultdict(int)
seen_IPR_set = set()

for line in interpro_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[0]
    IPR_set = set(split_line[1].split(','))
    if gene_id in set1_genes:
        print(f"Matched gene in Set1: {gene_id}")
        for IPR in IPR_set:
            set1_dict[IPR] += 1
        set1_count += 1
        [seen_IPR_set.add(x) for x in IPR_set]
    elif gene_id in set2_genes:
        print(f"Matched gene in Set2: {gene_id}")
        for IPR in IPR_set:
            set2_dict[IPR] += 1
        set2_count += 1
        [seen_IPR_set.add(x) for x in IPR_set]
    else:
        print(f"No match for gene: {gene_id}")

print("Final counts:")
print(f"{Set1_name}: {set1_count}")
print(f"{Set2_name}: {set2_count}")


print(set1_count)
print(set2_count)

for IPR in seen_IPR_set:
    set1_IPR = set1_dict[IPR]
    set2_IPR = set2_dict[IPR]
    set1_others = set1_count - set1_IPR
    set2_others = set2_count - set2_IPR
    outline1 = "\t".join([str(IPR), str(set1_IPR), str(set2_IPR)]) + "\n"
    outline2 = "\t".join(["Other genes", str(set1_others), str(set2_others)]) + "\n"

    outfile = "".join([conf.outdir, "/", IPR, "_fischertable.txt"])
    with open(outfile, 'w') as o:
        o.write("".join([outline1, outline2]))
