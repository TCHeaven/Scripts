#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob
import os
import argparse
import matplotlib.pyplot as plt

from sklearn.metrics import pairwise_distances
from sklearn.decomposition import PCA

# -------------------------
# Arguments
# -------------------------
parser = argparse.ArgumentParser(description="PCoA (EMU + NanoVI) with metadata colors")

parser.add_argument("-i", "--input", required=True, help="Input directory with *_rel-abundance.tsv files")
parser.add_argument("-o", "--output", required=True, help="Output directory")

parser.add_argument("-l", "--level", default="genus", help="Taxonomic level (genus/species/family etc.)")

parser.add_argument("--metadata", required=False, help="CSV with sample,group,color columns")

parser.add_argument("--sample_order", nargs="+", help="Optional sample order")
parser.add_argument("--sample_order_file", help="File with sample order (one per line)")

args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)

# -------------------------
# Load abundance tables
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
# Build sample × taxa matrix
# -------------------------
agg = data.groupby(["sample", args.level])[count_col].sum().reset_index()

matrix = agg.pivot(index="sample", columns=args.level, values=count_col).fillna(0)

# -------------------------
# Apply sample order (optional)
# -------------------------
sample_order = None

if args.sample_order_file:
    with open(args.sample_order_file) as f:
        sample_order = [x.strip() for x in f if x.strip()]
elif args.sample_order:
    sample_order = args.sample_order

if sample_order:
    sample_order = [s for s in sample_order if s in matrix.index]
    remaining = [s for s in matrix.index if s not in sample_order]
    matrix = matrix.loc[sample_order + remaining]

# -------------------------
# Normalize to relative abundance
# -------------------------
matrix = matrix.div(matrix.sum(axis=1), axis=0).fillna(0)

# -------------------------
# Bray-Curtis distance
# -------------------------
dist = pairwise_distances(matrix, metric="braycurtis")

# -------------------------
# PCoA (PCA on distance matrix approximation)
# -------------------------
coords = PCA(n_components=2).fit_transform(dist)

# -------------------------
# Load metadata (optional)
# -------------------------
groups = None
colors = None

if args.metadata:
    meta = pd.read_csv(args.metadata)

    if not {"sample", "group", "color"}.issubset(meta.columns):
        raise ValueError("Metadata must contain: sample, group, color columns")

    meta = meta.set_index("sample")
    meta = meta.loc[matrix.index]

    groups = meta["group"].values
    colors = meta["color"].values

# -------------------------
# Plot
# -------------------------
plt.figure(figsize=(8, 6))

if colors is not None:
    # plot each sample with its own metadata color
    for i, sample in enumerate(matrix.index):
        plt.scatter(
            coords[i, 0],
            coords[i, 1],
            color=colors[i],
            s=80
        )

    # clean legend (unique groups only)
    for g in sorted(set(groups)):
        idx = np.where(groups == g)[0][0]
        plt.scatter([], [], color=colors[idx], label=g)

    plt.legend(title="Group")

else:
    plt.scatter(coords[:, 0], coords[:, 1], color="black", s=80)

# sample labels
for i, s in enumerate(matrix.index):
    plt.text(coords[i, 0], coords[i, 1], s, fontsize=8)

plt.xlabel("PCoA1")
plt.ylabel("PCoA2")
plt.title(f"PCoA (Bray-Curtis) - {args.level}")

plt.tight_layout()

# -------------------------
# Save outputs
# -------------------------
base = os.path.join(args.output, f"pcoa_{args.level}")

plt.savefig(base + ".png", dpi=300)
plt.savefig(base + ".svg", format="svg")

print(f"Saved: {base}.png and {base}.svg")