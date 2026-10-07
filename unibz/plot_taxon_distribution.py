#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt


def plot_taxon_distribution(input_file, taxon, output_file):
    """
    Read a tab-separated Kraken-style output file, extract values
    for a specified taxon, bin the values, and create a bar plot.
    """
    # Read TSV; no header assumed
    df = pd.read_csv(input_file, sep="\t", header=None)
    # Check that the expected columns exist
    if df.shape[1] < 4:
        raise ValueError(
            f"Input file has only {df.shape[1]} columns; "
            "at least 4 columns are required."
        )
    # Select the target taxon and convert column 4 to numeric
    values = pd.to_numeric(
        df.loc[df[2] == taxon, 3],
        errors="coerce"
    ).dropna()
    if values.empty:
        raise ValueError(
            f"No numeric values found for taxon: {taxon}"
        )
    # Define bins:
    # 0–1000 in steps of 10
    # 10000–50000 in steps of 1000
    bins = list(range(0, 1001, 10)) + list(range(10000, 50001, 1000))
    # Bin values
    binned = pd.cut(
        values,
        bins=bins,
        right=False,
        include_lowest=True
    )
    counts = binned.value_counts().sort_index()
    edges = counts.index
    # Plot
    plt.figure(figsize=(16, 6))
    plt.bar(
        range(len(counts)),
        counts.values,
        width=1
    )
    plt.xticks(
        range(len(counts)),
        [str(int(interval.left)) for interval in edges],
        rotation=90
    )
    plt.xlabel("Column 4 value")
    plt.ylabel("Number of reads")
    plt.title(taxon)
    plt.tight_layout()
    # Save plot
    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()
    print(f"Taxon: {taxon}")
    print(f"Input: {input_file}")
    print(f"Reads plotted: {len(values)}")
    print(f"Output: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot the distribution of values in column 4 "
            "for a specified taxon."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input tab-separated output_nt.txt file"
    )
    parser.add_argument(
        "-t",
        "--taxon",
        required=True,
        help="Taxon name exactly as it appears in column 3"
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output PNG filename"
    )
    args = parser.parse_args()
    plot_taxon_distribution(
        input_file=args.input,
        taxon=args.taxon,
        output_file=args.output
    )


if __name__ == "__main__":
    main()
