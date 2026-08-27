#!/usr/bin/env python3

import argparse
import pandas as pd
import os


parser = argparse.ArgumentParser(
    description="Collapse ODGI path haplotype matrix into sample node counts and presence/absence matrix"
)

parser.add_argument(
    "--input",
    required=True,
    help="Output from: odgi paths -H"
)

parser.add_argument(
    "--output",
    required=True,
    help="Output sample x node presence/absence matrix"
)

args = parser.parse_args()


# ---------------------------------------------------------
# Read ODGI haplotype matrix
# ---------------------------------------------------------

df = pd.read_csv(
    args.input,
    sep="\t"
)


# ---------------------------------------------------------
# Identify node columns
# ---------------------------------------------------------

node_cols = [
    c for c in df.columns
    if c.startswith("node.")
]

if len(node_cols) == 0:
    raise ValueError(
        "No node columns found. Check that input is from 'odgi paths -H'"
    )


# ---------------------------------------------------------
# Extract sample name from path name
#
# Example:
# AT1_AO_11_ET_fasta#0#AT1_AO_11_ET_fasta_1#0
#
# becomes:
# AT1_AO_11_ET_fasta
# ---------------------------------------------------------

df["sample"] = (
    df["path.name"]
    .str.split("#")
    .str[0]
)


# ---------------------------------------------------------
# Convert node values to numeric
# (ODGI haplotype matrix should already be 0/1,
# but this protects against formatting issues)
# ---------------------------------------------------------

df[node_cols] = df[node_cols].apply(
    pd.to_numeric,
    errors="coerce"
).fillna(0)


# ---------------------------------------------------------
# Collapse contigs within each sample
#
# SUM keeps information about multiple occurrences
#
# Example:
#
# contig1 node100 = 1
# contig2 node100 = 1
#
# sample node100 = 2
# ---------------------------------------------------------

collapsed_counts = (
    df
    .groupby("sample")[node_cols]
    .sum()
)


# ---------------------------------------------------------
# Write node occurrence counts
# ---------------------------------------------------------

base, ext = os.path.splitext(args.output)

count_output = base + "_counts.tsv"

collapsed_counts.to_csv(
    count_output,
    sep="\t"
)


# ---------------------------------------------------------
# Convert counts to presence / absence
#
# >0 = present
#  0 = absent
# ---------------------------------------------------------

collapsed_presence = (
    collapsed_counts > 0
).astype(int)


# ---------------------------------------------------------
# Write presence/absence matrix
# ---------------------------------------------------------

collapsed_presence.to_csv(
    args.output,
    sep="\t"
)


# ---------------------------------------------------------
# Summary
# ---------------------------------------------------------

print(
    f"Input paths: {len(df)}"
)

print(
    f"Collapsed samples: {len(collapsed_counts)}"
)

print(
    f"Graph nodes: {len(node_cols)}"
)

print(
    f"Count matrix written: {count_output}"
)

print(
    f"Presence/absence matrix written: {args.output}"
)


# Show node counts per sample
print("\nTotal nodes per sample:")
print(
    collapsed_presence.sum(axis=1)
)


# Show duplicated nodes
duplicates = (
    collapsed_counts > 1
).sum(axis=1)

print("\nSamples containing duplicated nodes:")
print(duplicates)

