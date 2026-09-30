#!/usr/bin/env python3

"""
Plot repeat coverage across a genome from a RepeatMasker GFF3 file.

Example:
    python repeat_density.py \
        -gff GCF_000026205.1_Phytoplasma_mali.fasta.out.gff \
        -o P_mali_repeat_density \
        --window 1000

Outputs:
    <prefix>.png
    <prefix>.pdf
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot RepeatMasker repeat coverage across a genome."
    )
    parser.add_argument(
        "-gff",
        "--gff",
        required=True,
        help="Input RepeatMasker GFF3 file"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="repeat_density",
        help="Output filename prefix [default: repeat_density]"
    )
    parser.add_argument(
        "--window",
        type=int,
        default=1000,
        help="Window size in bp [default: 1000]"
    )
    parser.add_argument(
        "--genome-length",
        type=int,
        default=None,
        help="Genome length in bp. If omitted, read from ##sequence-region."
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="PNG resolution [default: 300]"
    )
    return parser.parse_args()


def get_genome_length(gff, supplied_length=None):
    if supplied_length is not None:
        return supplied_length
    with open(gff) as f:
        for line in f:
            if line.startswith("##sequence-region"):
                fields = line.strip().split()
                # GFF3:
                # ##sequence-region seqid start end
                if len(fields) >= 4:
                    start = int(fields[2])
                    end = int(fields[3])
                    return end - start + 1
    raise ValueError(
        "Could not determine genome length from GFF. "
        "Use --genome-length."
    )


def read_gff(gff):
    df = pd.read_csv(
        gff,
        sep="\t",
        comment="#",
        header=None
    )
    if df.empty:
        raise ValueError("No GFF features found.")
    if df.shape[1] != 9:
        raise ValueError(
            f"Expected 9 GFF columns, found {df.shape[1]}."
        )
    df.columns = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    ]
    # Keep RepeatMasker repeat annotations
    df = df[df["type"] == "dispersed_repeat"].copy()
    if df.empty:
        raise ValueError(
            "No 'dispersed_repeat' features found in GFF."
        )
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def calculate_repeat_coverage(df, genome_length, window):
    # Per-base Boolean array.
    #
    # This prevents overlapping RepeatMasker annotations
    # from being counted multiple times.
    repeat_mask = np.zeros(genome_length, dtype=bool)
    for start, end in zip(df["start"], df["end"]):
        # GFF coordinates are 1-based inclusive.
        # Python arrays are 0-based.
        start = max(1, start)
        end = min(genome_length, end)
        repeat_mask[start - 1:end] = True
    n_windows = int(np.ceil(genome_length / window))
    repeat_bp = []
    window_lengths = []
    for i in range(n_windows):
        start = i * window
        end = min((i + 1) * window, genome_length)
        repeat_bp.append(
            repeat_mask[start:end].sum()
        )
        window_lengths.append(
            end - start
        )
    repeat_bp = np.array(repeat_bp)
    window_lengths = np.array(window_lengths)
    repeat_percent = (
        repeat_bp / window_lengths * 100
    )
    midpoints = (
        np.arange(n_windows) * window
        + window_lengths / 2
    )
    return midpoints, repeat_bp, repeat_percent, repeat_mask


def make_plot(
    midpoints,
    repeat_percent,
    genome_length,
    window,
    output,
    dpi
):
    fig, ax = plt.subplots(
        figsize=(15, 5)
    )
    ax.bar(
        midpoints / 1000,
        repeat_percent,
        width=window / 1000,
        align="center"
    )
    ax.set_xlim(
        0,
        genome_length / 1000
    )
    ax.set_xlabel(
        "Genome position (kb)"
    )
    ax.set_ylabel(
        "Repeat coverage (%)"
    )
    ax.set_title(
        "Repeat distribution across genome"
    )
    # Major ticks every 50 kb
    tick_step = 50
    ticks = np.arange(
        0,
        genome_length / 1000 + tick_step,
        tick_step
    )
    ax.set_xticks(ticks)
    ax.grid(
        axis="y",
        alpha=0.3
    )
    plt.tight_layout()
    # PNG
    png_file = output + ".png"
    plt.savefig(
        png_file,
        dpi=dpi,
        bbox_inches="tight"
    )
    # PDF
    pdf_file = output + ".pdf"
    plt.savefig(
        pdf_file,
        bbox_inches="tight"
    )
    plt.close()
    return png_file, pdf_file


def main():
    args = parse_args()
    # -----------------------------------------------------
    # Check input
    # -----------------------------------------------------
    if not os.path.isfile(args.gff):
        raise FileNotFoundError(
            f"Input GFF not found: {args.gff}"
        )
    if args.window <= 0:
        raise ValueError(
            "Window size must be greater than zero."
        )
    # -----------------------------------------------------
    # Genome length
    # -----------------------------------------------------
    genome_length = get_genome_length(
        args.gff,
        args.genome_length
    )
    print(
        f"Genome length: {genome_length:,} bp"
    )
    # -----------------------------------------------------
    # Read GFF
    # -----------------------------------------------------
    df = read_gff(args.gff)
    print(
        f"Repeat features: {len(df):,}"
    )
    # -----------------------------------------------------
    # Calculate coverage
    # -----------------------------------------------------
    (
        midpoints,
        repeat_bp,
        repeat_percent,
        repeat_mask
    ) = calculate_repeat_coverage(
        df,
        genome_length,
        args.window
    )
    total_repeat_bp = repeat_mask.sum()
    total_repeat_percent = (
        total_repeat_bp /
        genome_length *
        100
    )
    print(
        f"Unique repeat-covered bases: "
        f"{total_repeat_bp:,}"
    )
    print(
        f"Genome repeat coverage: "
        f"{total_repeat_percent:.2f}%"
    )
    # -----------------------------------------------------
    # Plot
    # -----------------------------------------------------
    png_file, pdf_file = make_plot(
        midpoints,
        repeat_percent,
        genome_length,
        args.window,
        args.output,
        args.dpi
    )
    print()
    print("Output:")
    print(f"  {png_file}")
    print(f"  {pdf_file}")


if __name__ == "__main__":
    main()
