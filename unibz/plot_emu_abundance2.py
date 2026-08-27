#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import glob
import os
import numpy as np

# -------------------------
# Argument parsing
# -------------------------
parser = argparse.ArgumentParser(description="Plot EMU relative abundance as stacked barplot")

parser.add_argument(
    "-i", "--input",
    required=True,
    help="Directory containing *_rel-abundance.tsv files"
)

parser.add_argument(
    "-l", "--level",
    default="genus",
    help="Taxonomic level (e.g., species, genus, family)"
)

parser.add_argument(
    "-n", "--top_n",
    type=int,
    default=15,
    help="Number of top taxa to plot (others grouped as 'Other')"
)

parser.add_argument(
    "-r", "--order",
    nargs="+",
    help="Optional sample order (space-separated list)"
)

parser.add_argument(
    "--output",
    default="emu_stacked_barplot.svg",
    help="Output plot filename (SVG or PNG)"
)

args = parser.parse_args()

# -------------------------
# Load files
# -------------------------
files = glob.glob(os.path.join(args.input, "*_rel-abundance.tsv"))
if not files:
    raise ValueError("No *_rel-abundance.tsv files found in the input directory.")

dfs = []
for f in files:
    df = pd.read_csv(f, sep="\t")
    sample = os.path.basename(f).replace("_rel-abundance.tsv", "")
    df["sample"] = sample
    dfs.append(df)

data = pd.concat(dfs, ignore_index=True)

# -------------------------
# Check taxonomic level
# -------------------------
if args.level not in data.columns:
    raise ValueError(f"Taxonomic level '{args.level}' not found in EMU output columns.")

# -------------------------
# Aggregate at chosen level
# -------------------------
agg = data.groupby(["sample", args.level])["abundance"].sum().reset_index()

# -------------------------
# Select top N taxa and group others as "Other"
# -------------------------
top_taxa = agg.groupby(args.level)["abundance"].sum().nlargest(args.top_n).index
agg[args.level] = agg[args.level].where(agg[args.level].isin(top_taxa), "Other")

# Re-aggregate after grouping "Other"
agg = agg.groupby(["sample", args.level])["abundance"].sum().reset_index()

# -------------------------
# Pivot table
# -------------------------
pivot = agg.pivot_table(
    index="sample",
    columns=args.level,
    values="abundance",
    fill_value=0
)

# -------------------------
# Sample order
# -------------------------
if args.order:
    pivot = pivot.reindex(args.order)
else:
    pivot = pivot.sort_index()

# -------------------------
# Stack order
# -------------------------
taxa_order = pivot.sum().sort_values(ascending=False).index.tolist()

if "Other" in taxa_order:
    taxa_order.remove("Other")

taxa_order = taxa_order[::-1]

if "Other" in pivot.columns:
    taxa_order = ["Other"] + taxa_order

pivot = pivot[taxa_order]

# -------------------------
# Normalize to relative abundance
# -------------------------
pivot = pivot.div(pivot.sum(axis=1), axis=0)

# -------------------------
# COLOR SYSTEM (FINAL FIX)
# -------------------------
n_taxa = len(pivot.columns)

# 1. EXACT first 9 matplotlib default colors
base_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
first_colors = base_colors[:9]

# 2. Grey for "Other"
other_color = "#BDBDBD"

# 3. Remaining taxa (everything after first 9 excluding "Other")
taxa = list(pivot.columns)

# assign first 9 in order of appearance
color_map = {}

i_default = 0
remaining_taxa = []

for t in taxa:
    if t == "Other":
        continue
    if i_default < 9:
        color_map[t] = first_colors[i_default]
        i_default += 1
    else:
        remaining_taxa.append(t)

# 4. Unique HSV colors for remaining taxa (NO repetition possible)
n_remaining = len(remaining_taxa)

if n_remaining > 0:
    hsv_colors = plt.cm.hsv(np.linspace(0, 1, n_remaining, endpoint=False))
    for t, c in zip(remaining_taxa, hsv_colors):
        color_map[t] = c

# 5. Force Other = grey
if "Other" in pivot.columns:
    color_map["Other"] = other_color

# -------------------------
# PLOT
# -------------------------
fig, ax = plt.subplots(figsize=(16, 6))

pivot.plot(
    kind="bar",
    stacked=True,
    width=0.9,
    ax=ax,
    color=[color_map[col] for col in pivot.columns]
)

plt.ylabel("Relative abundance")
plt.xlabel("Sample")
plt.xticks(rotation=45, ha="right")

# Legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles[::-1],
    labels[::-1],
    title=args.level,
    bbox_to_anchor=(1.05, 1),
    loc="upper left"
)

plt.tight_layout()

# -------------------------
# Save
# -------------------------
outdir = os.path.dirname(args.output)
if outdir:
    os.makedirs(outdir, exist_ok=True)

plt.savefig(args.output, format="svg")
print(f"Stacked bar plot saved to: {args.output}")
