#!/usr/bin/env python3

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse

# ----------------------------
# Arguments
# ----------------------------
parser = argparse.ArgumentParser(description="Stacked barplot for taxa abundance")

parser.add_argument("-i", "--input", required=True, help="Input TSV file")
parser.add_argument("-o", "--output", default="stacked_barplot.svg", help="Output plot")
parser.add_argument("-n", "--top_n", type=int, default=10, help="Top N taxa")
parser.add_argument("--sample_order", nargs="+", help="Custom sample order")
parser.add_argument("--tax_level", default="Genus",
                    choices=["Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"])

args = parser.parse_args()

# ----------------------------
# Load data
# ----------------------------
df = pd.read_csv(args.input, sep="\t")

# ----------------------------
# Remove 'total' column if present
# ----------------------------
if "total" in df.columns:
    df = df.drop(columns=["total"])

# ----------------------------
# Split taxonomy column
# ----------------------------
tax_split = df["tax"].str.split(";", expand=True)

all_levels = ["Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"]

n_levels = tax_split.shape[1]

if n_levels > len(all_levels):
    raise ValueError(f"Too many taxonomy levels ({n_levels}) in input")

levels_used = all_levels[:n_levels]
tax_split.columns = levels_used

# Merge back into dataframe
df = pd.concat([df, tax_split], axis=1)

# ----------------------------
# Check requested tax level exists
# ----------------------------
if args.tax_level not in levels_used:
    raise ValueError(
        f"Tax level '{args.tax_level}' not found. Available: {levels_used}"
    )

# ----------------------------
# Select taxonomic level
# ----------------------------
df["Taxon"] = df[args.tax_level]

# Clean labels
df["Taxon"] = df["Taxon"].replace(["", "Unknown"], "Unclassified")

# ----------------------------
# Drop taxonomy columns
# ----------------------------
df = df.drop(columns=["tax"] + levels_used)

# ----------------------------
# Set index
# ----------------------------
df = df.set_index("Taxon")

# ----------------------------
# Convert to long format
# ----------------------------
df_long = df.reset_index().melt(
    id_vars="Taxon",
    var_name="Sample",
    value_name="Abundance"
)

# ----------------------------
# Aggregate counts
# ----------------------------
df_long = df_long.groupby(["Sample", "Taxon"], as_index=False).sum()

# ----------------------------
# Convert to relative abundance
# ----------------------------
df_long["Abundance"] = df_long.groupby("Sample")["Abundance"].transform(
    lambda x: x / x.sum()
)

# ----------------------------
# Select top N taxa (global)
# ----------------------------
top_taxa = (
    df_long.groupby("Taxon")["Abundance"]
    .sum()
    .sort_values(ascending=False)
    .head(args.top_n)
    .index
)

# Group others
df_long["Taxon"] = df_long["Taxon"].where(df_long["Taxon"].isin(top_taxa), "Other")

# Re-aggregate after grouping
df_long = df_long.groupby(["Sample", "Taxon"], as_index=False).sum()

# ----------------------------
# Pivot for plotting
# ----------------------------
df_pivot = df_long.pivot(index="Sample", columns="Taxon", values="Abundance").fillna(0)

# Convert to %
df_pivot = df_pivot * 100

# ----------------------------
# Apply sample order if provided
# ----------------------------
if args.sample_order:
    df_pivot = df_pivot.loc[args.sample_order]

# ----------------------------
# Order taxa (Other always at bottom)
# ----------------------------

# compute global abundance per taxon
taxa_order = (
    df_long.groupby("Taxon")["Abundance"]
    .sum()
    .sort_values(ascending=True)   # <-- IMPORTANT: ascending
    .index
    .tolist()
)

# remove "Other" from ranking
taxa_order = [t for t in taxa_order if t != "Other"]

# build final column order
cols = []

# bottom of stack
if "Other" in df_pivot.columns:
    cols.append("Other")

# then least → most abundant (so most ends up on top)
cols.extend([t for t in taxa_order if t in df_pivot.columns])

df_pivot = df_pivot[cols]

# ----------------------------
# Save stack order for legend sync
# ----------------------------

stack_order = df_pivot.columns.tolist()

# ----------------------------
# Plot
# ----------------------------
ax = df_pivot.plot(kind="bar", stacked=True, figsize=(16, 6))

ax.set_xlabel("Sample")
ax.set_ylabel("Relative abundance (%)")
ax.set_title(f"Top {args.top_n} taxa at {args.tax_level} level")

plt.xticks(rotation=45, ha="right")

# ----------------------------
# Move legend outside
# ----------------------------
handles, labels = ax.get_legend_handles_labels()

handle_map = dict(zip(labels, handles))

legend_order = stack_order[::-1]

ordered_handles = []
ordered_labels = []

for taxon in legend_order:
    if taxon in handle_map:
        ordered_handles.append(handle_map[taxon])
        ordered_labels.append(taxon)

ax.legend(
    ordered_handles,
    ordered_labels,
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False
)

plt.tight_layout()

plt.savefig(args.output)
