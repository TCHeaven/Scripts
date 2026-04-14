#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import argparse

# -------------------------
# Arguments
# -------------------------
parser = argparse.ArgumentParser(description="Reusable microbiome heatmap (EMU + NanoVI)")

parser.add_argument("-i", "--input", required=True, help="Input directory with *_rel-abundance.tsv files")
parser.add_argument("-o", "--output", required=True, help="Output directory")

parser.add_argument("-l", "--level", default="genus", help="Taxonomic level (species, genus, family, etc.)")

parser.add_argument("--log", action="store_true", help="Apply log1p transform")

parser.add_argument("--top_n", type=int, default=30, help="Keep top N taxa by abundance")

parser.add_argument("--sample_order", nargs="+", help="Sample order (space-separated)")

parser.add_argument("--sample_order_file", help="File with sample order (one per line)")

args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)

# -------------------------
# Load all files
# -------------------------
files = glob.glob(os.path.join(args.input, "*_rel-abundance.tsv"))

if not files:
    raise RuntimeError("No *_rel-abundance.tsv files found")

dfs = []

for f in files:
    df = pd.read_csv(f, sep="\t")
    sample = os.path.basename(f).replace("_rel-abundance.tsv", "")
    df["sample"] = sample
    dfs.append(df)

data = pd.concat(dfs, ignore_index=True)

# -------------------------
# Detect abundance column
# -------------------------
if "estimated_counts" in data.columns:
    count_col = "estimated_counts"
elif "estimated counts" in data.columns:
    count_col = "estimated counts"
else:
    raise RuntimeError("No abundance column found")

# -------------------------
# Ensure taxonomy column exists
# -------------------------
if args.level not in data.columns:
    data[args.level] = "Unknown"

# -------------------------
# Aggregate
# -------------------------
agg = data.groupby(["sample", args.level])[count_col].sum().reset_index()

matrix = agg.pivot(index=args.level, columns="sample", values=count_col).fillna(0)

# -------------------------
# Apply sample order (NEW)
# -------------------------
sample_order = None

if args.sample_order_file:
    with open(args.sample_order_file) as f:
        sample_order = [x.strip() for x in f if x.strip()]

elif args.sample_order:
    sample_order = args.sample_order

if sample_order:
    sample_order = [s for s in sample_order if s in matrix.columns]
    remaining = [s for s in matrix.columns if s not in sample_order]
    matrix = matrix[sample_order + remaining]

# -------------------------
# Filter top taxa
# -------------------------
top_taxa = matrix.sum(axis=1).sort_values(ascending=False).head(args.top_n).index
matrix = matrix.loc[top_taxa]

# -------------------------
# Transform
# -------------------------
if args.log:
    matrix = np.log1p(matrix)

# -------------------------
# Plot
# -------------------------
plt.figure(figsize=(12, max(6, len(matrix) * 0.3)))

sns.heatmap(matrix, cmap="viridis")

plt.title(f"Microbiome Heatmap ({args.level})")
plt.xlabel("Samples")
plt.ylabel(args.level)

plt.xticks(rotation=90)

plt.tight_layout()

# -------------------------
# Save outputs
# -------------------------
base = os.path.join(args.output, f"heatmap_{args.level}")

plt.savefig(base + ".png", dpi=300)
plt.savefig(base + ".svg", format="svg")

print(f"Saved: {base}.png and {base}.svg")