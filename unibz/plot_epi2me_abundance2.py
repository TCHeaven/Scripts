#!/usr/bin/env python3

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse

# ----------------------------
# Arguments
# ----------------------------
parser = argparse.ArgumentParser(description="Stacked barplot for taxa abundance (EMU-style)")

parser.add_argument("-i", "--input", required=True, help="Input TSV file")
parser.add_argument("-o", "--output", default="stacked_barplot.svg", help="Output plot")
parser.add_argument("-n", "--top_n", type=int, default=10, help="Top N taxa")
parser.add_argument("--sample_order", nargs="+", help="Custom sample order")

parser.add_argument(
    "--tax_level",
    default="Genus",
    choices=["Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"]
)

args = parser.parse_args()

# ----------------------------
# Load data
# ----------------------------
df = pd.read_csv(args.input, sep="\t")

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

df = pd.concat([df, tax_split], axis=1)

if args.tax_level not in levels_used:
    raise ValueError(
        f"Tax level '{args.tax_level}' not found. Available: {levels_used}"
    )

df["Taxon"] = df[args.tax_level]
df["Taxon"] = df["Taxon"].replace(["", "Unknown"], "Unclassified")

df = df.drop(columns=["tax"] + levels_used)

# ----------------------------
# Long format
# ----------------------------
df = df.set_index("Taxon")

df_long = df.reset_index().melt(
    id_vars="Taxon",
    var_name="Sample",
    value_name="Abundance"
)

df_long = df_long.groupby(["Sample", "Taxon"], as_index=False).sum()

# ----------------------------
# Top N taxa + Other grouping
# ----------------------------
top_taxa = (
    df_long.groupby("Taxon")["Abundance"]
    .sum()
    .sort_values(ascending=False)
    .head(args.top_n)
    .index
)

df_long["Taxon"] = df_long["Taxon"].where(df_long["Taxon"].isin(top_taxa), "Other")

df_long = df_long.groupby(["Sample", "Taxon"], as_index=False).sum()

# ----------------------------
# Pivot (EMU style)
# ----------------------------
df_pivot = df_long.pivot(
    index="Sample",
    columns="Taxon",
    values="Abundance"
).fillna(0)

# ----------------------------
# Sample order
# ----------------------------
if args.sample_order:
    df_pivot = df_pivot.reindex(args.sample_order)

# ----------------------------
# EMU normalization (IMPORTANT)
# ----------------------------
df_pivot = df_pivot.div(df_pivot.sum(axis=1), axis=0)

# ----------------------------
# Taxa ordering (EMU style)
# ----------------------------
taxa_order = df_pivot.sum().sort_values(ascending=False).index.tolist()

# remove Other from ranking
if "Other" in taxa_order:
    taxa_order.remove("Other")

# least → most (so most ends on top of stack)
taxa_order = taxa_order[::-1]

# ensure Other is at bottom
if "Other" in df_pivot.columns:
    taxa_order = ["Other"] + taxa_order

df_pivot = df_pivot[taxa_order]

# ----------------------------
# Plot
# ----------------------------
fig, ax = plt.subplots(figsize=(16, 6))

df_pivot.plot(
    kind="bar",
    stacked=True,
    width=0.9,
    ax=ax
)

ax.set_xlabel("Sample")
ax.set_ylabel("Relative abundance")
ax.set_title(f"Top {args.top_n} taxa at {args.tax_level} level")

plt.xticks(rotation=45, ha="right")

# ----------------------------
# Legend (EMU style reversed)
# ----------------------------
handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles[::-1],
    labels[::-1],
    title=args.tax_level,
    bbox_to_anchor=(1.05, 1),
    loc="upper left"
)

# ----------------------------
# Layout + save
# ----------------------------
plt.tight_layout()
plt.savefig(args.output, format="svg", bbox_inches="tight")
