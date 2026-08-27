#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt

from upsetplot import UpSet, from_indicators


parser = argparse.ArgumentParser(
    description="Create UpSet plot from Orthofinder orthogroup presence matrix"
)

parser.add_argument(
    "--input",
    required=True,
    help="Orthofinder orthogroup_presence.tsv"
)

parser.add_argument(
    "--output",
    required=True,
    help="Output PNG file"
)

args = parser.parse_args()


print(f"Reading {args.input}")

# Read table
df = pd.read_csv(
    args.input,
    sep="\t",
    dtype=str
)


# Remove orthogroup ID
orthogroups = df["Orthogroup"]

presence = df.drop(columns=["Orthogroup"])


# Convert True/False strings to boolean
presence = presence.applymap(
    lambda x: True if str(x).strip().lower() == "true" else False
)


print("\nOrthogroups present per genome:")
print(presence.sum())

print("\nTotal orthogroups:")
print(len(presence))


# Convert to UpSet format
upset_data = from_indicators(
    presence.columns,
    presence
)


plt.figure(figsize=(12,8))


UpSet(
    upset_data,
    show_counts=True,
    sort_by="cardinality"
).plot()


plt.savefig(
    args.output,
    dpi=300,
    bbox_inches="tight"
)

plt.close()


print(f"\nWritten: {args.output}")
