#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description="Plot PCA results from bcftools +pca output"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="PCA output file from bcftools +pca"
    )
    parser.add_argument(
        "-o", "--output",
        default="pca_plot.png",
        help="Output image filename"
    )
    parser.add_argument(
        "--pcx",
        default="PC1",
        help="PC for x-axis (default: PC1)"
    )
    parser.add_argument(
        "--pcy",
        default="PC2",
        help="PC for y-axis (default: PC2)"
    )
    parser.add_argument(
        "--labels",
        action="store_true",
        help="Label samples on plot"
    )
    parser.add_argument(
        "--fontsize",
        type=int,
        default=8,
        help="Sample label font size"
    )
    parser.add_argument(
        "--width",
        type=float,
        default=8,
        help="Figure width"
    )
    parser.add_argument(
        "--height",
        type=float,
        default=6,
        help="Figure height"
    )
    args = parser.parse_args()
    # Read bcftools PCA output
    df = pd.read_csv(
        args.input,
        sep="\t",
        comment="#",
        header=None
    )
    # Assign column names
    n_pc = df.shape[1] - 1
    columns = ["sample"] + [
        f"PC{i}" for i in range(1, n_pc + 1)
    ]
    df.columns = columns
    # Check requested PCs exist
    for pc in [args.pcx, args.pcy]:
        if pc not in df.columns:
            raise ValueError(
                f"{pc} not found. Available PCs: {', '.join(df.columns[1:])}"
            )
    # Plot
    plt.figure(
        figsize=(args.width, args.height)
    )
    plt.scatter(
        df[args.pcx],
        df[args.pcy]
    )
    if args.labels:
        for _, row in df.iterrows():
            plt.text(
                row[args.pcx],
                row[args.pcy],
                row["sample"],
                fontsize=args.fontsize
            )
    plt.xlabel(args.pcx)
    plt.ylabel(args.pcy)
    plt.title(
        f"PCA: {args.pcx} vs {args.pcy}"
    )
    plt.tight_layout()
    plt.savefig(
        args.output,
        dpi=300
    )
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
