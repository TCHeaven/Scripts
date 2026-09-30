#!/usr/bin/env python3

import argparse
import subprocess
import sys

import numpy as np
import matplotlib.pyplot as plt


def get_contig_length(bam, contig):
    """Get contig length from BAM header."""
    cmd = ["samtools", "view", "-H", bam]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True
    )
    for line in result.stdout.splitlines():
        if line.startswith("@SQ"):
            fields = dict(
                field.split(":", 1)
                for field in line.split("\t")[1:]
                if ":" in field
            )
            if fields.get("SN") == contig:
                return int(fields["LN"])
    raise ValueError(f"Contig '{contig}' not found in BAM header.")


def get_coverage(bam, contig, length):
    """
    Get per-base coverage using samtools depth.
    -a ensures positions with zero coverage are included.
    """
    cmd = [
        "samtools",
        "depth",
        "-a",
        "-r", contig,
        bam
    ]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True
    )
    coverage = np.zeros(length, dtype=np.int32)
    for line in result.stdout.splitlines():
        if not line:
            continue
        chrom, pos, depth = line.split("\t")
        pos = int(pos)
        # Convert 1-based SAM position to 0-based Python index
        coverage[pos - 1] = int(depth)
    return coverage


def find_intervals(mask):
    """Convert boolean mask into contiguous 1-based intervals."""
    intervals = []
    in_interval = False
    start = None
    for i, value in enumerate(mask):
        if value and not in_interval:
            start = i + 1
            in_interval = True
        elif not value and in_interval:
            end = i
            intervals.append((start, end))
            in_interval = False
    if in_interval:
        intervals.append((start, len(mask)))
    return intervals


def main():
    parser = argparse.ArgumentParser(
        description="Plot low-coverage positions from a BAM alignment."
    )
    parser.add_argument(
        "-b", "--bam",
        required=True,
        help="Input sorted/indexed BAM file"
    )
    parser.add_argument(
        "-c", "--contig",
        required=True,
        help="Reference contig/chromosome name"
    )
    parser.add_argument(
        "-t", "--threshold",
        type=int,
        default=5,
        help="Flag positions with coverage below this value (default: 5)"
    )
    parser.add_argument(
        "-o", "--output",
        default="coverage.png",
        help="Output plot filename (default: coverage.png)"
    )
    parser.add_argument(
        "--intervals",
        default=None,
        help="Optional output TSV for low-coverage intervals"
    )
    args = parser.parse_args()
    print(f"BAM:       {args.bam}")
    print(f"Contig:    {args.contig}")
    print(f"Threshold: coverage < {args.threshold}")
    length = get_contig_length(args.bam, args.contig)
    print(f"Length:    {length:,} bp")
    coverage = get_coverage(
        args.bam,
        args.contig,
        length
    )
    positions = np.arange(1, length + 1)
    low = coverage < args.threshold
    print(f"Low-coverage positions: {low.sum():,}")
    print(
        f"Fraction low coverage: "
        f"{100 * low.mean():.3f}%"
    )
    intervals = find_intervals(low)
    print(f"Low-coverage intervals: {len(intervals):,}")
    if intervals:
        print("\nFirst low-coverage intervals:")
        for start, end in intervals[:20]:
            print(
                f"  {start:,}\t{end:,}\t"
                f"{end - start + 1:,} bp"
            )
    # Write intervals if requested
    if args.intervals:
        with open(args.intervals, "w") as out:
            out.write("contig\tstart\tend\tlength\n")
            for start, end in intervals:
                out.write(
                    f"{args.contig}\t"
                    f"{start}\t"
                    f"{end}\t"
                    f"{end - start + 1}\n"
                )
        print(f"\nIntervals written to: {args.intervals}")
    # -----------------------------
    # Plot
    # -----------------------------
    fig, ax = plt.subplots(
        figsize=(16, 5)
    )
    ax.plot(
        positions / 1000,
        coverage,
        linewidth=0.5
    )
    # Highlight low-coverage positions
    ax.scatter(
        positions[low] / 1000,
        coverage[low],
        s=2
    )
    ax.axhline(
        args.threshold,
        linestyle="--",
        linewidth=1,
        label=f"threshold = {args.threshold}"
    )
    ax.set_xlabel("Genome position (kb)")
    ax.set_ylabel("Read depth")
    ax.set_title(
        f"Coverage across {args.contig}"
    )
    ax.set_xlim(
        0,
        length / 1000
    )
    ax.legend()
    plt.tight_layout()
    plt.savefig(
        args.output,
        dpi=300
    )
    print(f"\nPlot written to: {args.output}")


if __name__ == "__main__":
    main()
