#!/usr/bin/env python3
"""
make_qiime2_inputs.py

Convert a "samples x sequences" count table plus an ASV ID <-> sequence map
into QIIME 2-ready inputs:
  - ASV_table.tsv  (FeatureID x SampleID, counts)
  - ASVs.fasta     (rep sequences with FeatureID headers)

Example:
  python make_qiime2_inputs.py \
    --counts ASV_18S_asv.tsv \
    --map ASV_id_map.csv \
    --out-prefix run1 \
    --out-dir out_qiime2 \
    --merge-duplicates sum
"""

from __future__ import annotations
import argparse
import os
import sys
from typing import Dict, List, Tuple

import pandas as pd


def norm_seq(s: str) -> str:
    """Normalize a sequence string for matching."""
    return str(s).upper().strip().replace(" ", "")


def read_map(map_path: str, asv_col: str, seq_col: str) -> pd.DataFrame:
    m = pd.read_csv(map_path)
    if asv_col not in m.columns or seq_col not in m.columns:
        raise ValueError(
            f"Mapping file must contain columns '{asv_col}' and '{seq_col}'. "
            f"Found: {list(m.columns)}"
        )
    m = m[[asv_col, seq_col]].copy()
    m[asv_col] = m[asv_col].astype(str).str.strip()
    m[seq_col] = m[seq_col].astype(str).map(norm_seq)

    # Basic validation
    if m[asv_col].duplicated().any():
        dup = m.loc[m[asv_col].duplicated(), asv_col].iloc[0]
        raise ValueError(f"Duplicate ASV IDs in mapping file (e.g., '{dup}').")
    if m[seq_col].duplicated().any():
        # Multiple ASV IDs mapping to identical sequences can happen, but is unusual.
        # We'll keep the first and warn.
        pass

    return m


def read_counts(counts_path: str, sep: str) -> pd.DataFrame:
    # Read header as-is (sequence strings are column names)
    df = pd.read_csv(counts_path, sep=sep, dtype=str)
    if df.shape[1] < 2:
        raise ValueError("Counts table must have at least 2 columns (SampleID + features).")
    return df


def rename_seq_columns_to_asv(
    df: pd.DataFrame,
    seq_to_asv: Dict[str, str],
    sample_col: str | None,
) -> Tuple[pd.DataFrame, List[str], List[str]]:
    # Determine sample column
    if sample_col is None:
        sample_col = df.columns[0]
    if sample_col not in df.columns:
        raise ValueError(f"Sample column '{sample_col}' not found. Available: {list(df.columns)}")

    df = df.rename(columns={sample_col: "SampleID"}).copy()
    df["SampleID"] = df["SampleID"].astype(str).str.strip()

    # Map sequence columns -> ASV IDs
    missing: List[str] = []
    kept_cols: List[str] = ["SampleID"]
    asv_cols: List[str] = []

    for c in df.columns:
        if c == "SampleID":
            continue
        seq = norm_seq(c)
        asv = seq_to_asv.get(seq)
        if asv is None:
            missing.append(c)
            continue
        kept_cols.append(c)
        asv_cols.append(asv)

    # Keep only mapped columns
    df2 = df[kept_cols].copy()

    # Rename mapped columns to ASV IDs
    rename_map = {orig: asv for orig, asv in zip(kept_cols[1:], asv_cols)}
    df2 = df2.rename(columns=rename_map)

    return df2, missing, asv_cols


def merge_duplicate_asv_columns(df: pd.DataFrame, how: str) -> pd.DataFrame:
    """
    If multiple sequence-columns map to the same ASV ID (after renaming),
    merge them by sum/max/first.
    """
    cols = list(df.columns)
    if "SampleID" not in cols:
        raise ValueError("Expected 'SampleID' column after renaming.")

    feature_cols = [c for c in cols if c != "SampleID"]
    if len(feature_cols) == len(set(feature_cols)):
        return df  # no duplicates

    # Convert to numeric early
    for c in feature_cols:
        df[c] = pd.to_numeric(df[c], errors="raise").astype(int)

    if how == "sum":
        merged = df.groupby("SampleID", as_index=False).sum(numeric_only=True)
        # groupby doesn't merge duplicate *columns*; so we do it manually:
        merged = (
            df.set_index("SampleID")
              .groupby(level=0)
              .sum()
              .reset_index()
        )
        # Now collapse duplicate columns by grouping columns axis
        merged = (
            merged.set_index("SampleID")
                  .groupby(axis=1, level=0)
                  .sum()
                  .reset_index()
        )
        return merged

    if how == "max":
        merged = (
            df.set_index("SampleID")
              .groupby(axis=1, level=0)
              .max()
              .reset_index()
        )
        return merged

    if how == "first":
        # Keep the first occurrence of each duplicate column
        keep = ["SampleID"]
        seen = set()
        for c in feature_cols:
            if c not in seen:
                keep.append(c)
                seen.add(c)
        return df[keep].copy()

    raise ValueError(f"Unknown merge method: {how}")


def write_outputs(
    df_counts_samples_rows: pd.DataFrame,
    map_df: pd.DataFrame,
    asv_col: str,
    seq_col: str,
    out_dir: str,
    out_prefix: str,
) -> Tuple[str, str]:
    os.makedirs(out_dir, exist_ok=True)

    # Ensure numeric counts
    for c in df_counts_samples_rows.columns:
        if c == "SampleID":
            continue
        df_counts_samples_rows[c] = pd.to_numeric(df_counts_samples_rows[c], errors="raise").astype(int)

    # Transpose: rows=FeatureID (ASVs), cols=samples
    df_t = df_counts_samples_rows.set_index("SampleID").T
    df_t.index.name = "FeatureID"

    table_path = os.path.join(out_dir, f"{out_prefix}ASV_table.tsv")
    fasta_path = os.path.join(out_dir, f"{out_prefix}ASVs.fasta")

    df_t.to_csv(table_path, sep="\t")

    # Write FASTA in mapping order
    with open(fasta_path, "w") as f:
        for _, row in map_df.iterrows():
            asv = str(row[asv_col]).strip()
            seq = str(row[seq_col]).strip().upper()
            f.write(f">{asv}\n{seq}\n")

    return table_path, fasta_path


def main() -> int:
    p = argparse.ArgumentParser(
        description="Create QIIME2-ready ASV table TSV and FASTA from counts + ASV<->sequence map."
    )
    p.add_argument("--counts", required=True, help="Counts table (samples as rows, sequences as columns).")
    p.add_argument("--map", required=True, help="Mapping CSV with ASV IDs and sequences.")
    p.add_argument("--sep", default="\t", help="Counts table delimiter (default: tab).")
    p.add_argument("--sample-col", default=None, help="Name of the sample ID column (default: first column).")
    p.add_argument("--asv-col", default="ASV", help="Column name for ASV IDs in mapping file (default: ASV).")
    p.add_argument("--seq-col", default="Sequence", help="Column name for sequences in mapping file (default: Sequence).")

    p.add_argument(
        "--out-dir", default=".", help="Output directory (default: current directory)."
    )
    p.add_argument(
        "--out-prefix",
        default="",
        help="Prefix for output files (e.g. 'run1_' -> run1_ASV_table.tsv, run1_ASVs.fasta).",
    )

    p.add_argument(
        "--merge-duplicates",
        choices=["sum", "max", "first", "error"],
        default="error",
        help="How to handle duplicate ASV columns after mapping (default: error).",
    )

    args = p.parse_args()

    # Read inputs
    m = read_map(args.map, args.asv_col, args.seq_col)
    seq_to_asv = dict(zip(m[args.seq_col], m[args.asv_col]))

    df = read_counts(args.counts, args.sep)

    # Rename columns
    df2, missing, mapped_asvs = rename_seq_columns_to_asv(df, seq_to_asv, args.sample_col)

    # Report missing columns
    if missing:
        print(
            f"WARNING: {len(missing)} columns in counts table did not match any sequence in the mapping file.",
            file=sys.stderr,
        )
        print("First few unmapped column headers:", missing[:3], file=sys.stderr)

    # Handle duplicate columns after mapping
    feature_cols = [c for c in df2.columns if c != "SampleID"]
    if len(feature_cols) != len(set(feature_cols)):
        if args.merge_duplicates == "error":
            dups = pd.Series(feature_cols).value_counts()
            dups = dups[dups > 1].index.tolist()[:10]
            raise ValueError(
                "Duplicate ASV IDs detected after mapping (same ASV mapped from multiple sequence columns).\n"
                f"Examples: {dups}\n"
                "Re-run with --merge-duplicates sum|max|first to resolve automatically."
            )
        df2 = merge_duplicate_asv_columns(df2, args.merge_duplicates)

    # Write outputs
    table_path, fasta_path = write_outputs(
        df2, m, args.asv_col, args.seq_col, args.out_dir, args.out_prefix
    )

    print(f"Wrote: {table_path}")
    print(f"Wrote: {fasta_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
