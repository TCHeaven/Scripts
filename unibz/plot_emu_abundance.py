#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import glob
import os

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
pivot = agg.pivot_table(index="sample", columns=args.level, values="abundance", fill_value=0)

# -------------------------
# Sample order
# -------------------------
if args.order:
    pivot = pivot.reindex(args.order)
else:
    pivot = pivot.sort_index()

# -------------------------
# Stack order: Most abundant on top, Other at bottom
# -------------------------
taxa_order = pivot.sum().sort_values(ascending=False).index.tolist()
if "Other" in taxa_order:
    taxa_order.remove("Other")
taxa_order = taxa_order[::-1]  # least → most
if "Other" in pivot.columns:
    taxa_order = ["Other"] + taxa_order  # Other at bottom
pivot = pivot[taxa_order]

# -------------------------
# Normalize to relative abundance (0–1)
# -------------------------
pivot = pivot.div(pivot.sum(axis=1), axis=0)

# -------------------------
# Plot
# -------------------------
fig, ax = plt.subplots(figsize=(16, 6))  # wider figure for more samples
pivot.plot(kind="bar", stacked=True, width=0.9, ax=ax)

plt.ylabel("Relative abundance")
plt.xlabel("Sample")
plt.xticks(rotation=45, ha="right")

# Legend top → bottom
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], title=args.level,
          bbox_to_anchor=(1.05, 1), loc="upper left")

plt.tight_layout()

# Ensure output directory exists
outdir = os.path.dirname(args.output)
if outdir:
    os.makedirs(outdir, exist_ok=True)

# Save as SVG
plt.savefig(args.output, format="svg")
print(f"Stacked bar plot saved to: {args.output}")
