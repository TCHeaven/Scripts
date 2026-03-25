#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from typing import Dict, List, Tuple

ASV_RE = re.compile(r"^ASV(\d+)$")

def asv_sort_key(asv: str):
    m = ASV_RE.match(asv)
    return (0, int(m.group(1))) if m else (1, asv)

def parse_input(arg: str) -> Tuple[str, str, int]:
    # SOURCE:PATH:ASVCOL (ASVCOL is 1-based)
    try:
        source, path, col = arg.split(":", 2)
        return source, path, int(col) - 1
    except Exception:
        raise ValueError(f"Invalid --input '{arg}'. Use SOURCE:PATH:ASV_COLUMN")

def _fix_header_tokens(tokens: List[str]) -> List[str]:
    # Repair known BLAST header issue where evalue and bitscore get stuck together
    fixed: List[str] = []
    for t in tokens:
        if t == "evaluebitscore":
            fixed.extend(["evalue", "bitscore"])
        else:
            fixed.append(t)
    return fixed

def read_table(path: str, asv_col_index: int) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    """
    Reads TSV/CSV/whitespace tables with a header line.
    Key improvement: if whitespace-delimited, split with maxsplit so last column can contain spaces (e.g. BLAST stitle).
    Returns:
      - columns excluding ASV column
      - data: {ASV -> {col -> value}}
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    # Read header line raw
    with open(path, "r", newline="") as f:
        header_line = f.readline()
        if not header_line:
            raise ValueError(f"Empty file: {path}")

    header_line = header_line.rstrip("\n\r")

    # Choose delimiter mode
    # If any tabs exist, treat as TSV (safe for BLAST outfmt when tab-separated)
    # Else if .csv, comma
    # Else whitespace
    if "\t" in header_line:
        mode = "delim"
        delim = "\t"
        header = [h.strip() for h in header_line.split("\t")]
    elif path.lower().endswith(".csv"):
        mode = "csv"
        delim = ","
        # parse properly as CSV header
        with open(path, "r", newline="") as f:
            header = next(csv.reader(f, delimiter=delim))
        header = [h.strip() for h in header]
    else:
        mode = "ws"
        header = _fix_header_tokens(header_line.split())

    if asv_col_index >= len(header):
        raise ValueError(
            f"ASV column index {asv_col_index + 1} out of range for {path} (header has {len(header)} columns)"
        )

    out_cols = [h for i, h in enumerate(header) if i != asv_col_index]
    data: Dict[str, Dict[str, str]] = {}

    with open(path, "r", newline="") as f:
        # consume header line already read if using delim/ws; for csv we already parsed header separately
        if mode in ("delim", "ws"):
            _ = f.readline()

        if mode == "csv":
            reader = csv.reader(f, delimiter=delim)
            # header already consumed above; continue with remaining rows
            for row in reader:
                if not row or all(c.strip() == "" for c in row):
                    continue
                if len(row) < len(header):
                    row += [""] * (len(header) - len(row))

                asv = row[asv_col_index].strip()
                if not asv:
                    continue

                rec: Dict[str, str] = {}
                j = 0
                for i in range(len(header)):
                    if i == asv_col_index:
                        continue
                    rec[out_cols[j]] = row[i].strip() if i < len(row) else ""
                    j += 1
                data[asv] = rec

        elif mode == "delim":
            # delimiter split (TSV)
            for line in f:
                line = line.rstrip("\n\r")
                if not line:
                    continue
                row = line.split(delim)
                if len(row) < len(header):
                    row += [""] * (len(header) - len(row))
                # If row is longer than header, keep extras joined into the last column
                if len(row) > len(header):
                    row = row[: len(header) - 1] + [delim.join(row[len(header) - 1 :])]

                asv = row[asv_col_index].strip()
                if not asv:
                    continue

                rec: Dict[str, str] = {}
                j = 0
                for i in range(len(header)):
                    if i == asv_col_index:
                        continue
                    rec[out_cols[j]] = row[i].strip() if i < len(row) else ""
                    j += 1
                data[asv] = rec

        else:
            # whitespace split WITH maxsplit to preserve final free-text column (e.g. BLAST stitle)
            maxsplit = len(header) - 1
            for line in f:
                line = line.strip()
                if not line:
                    continue
                row = _fix_header_tokens(line.split(None, maxsplit=maxsplit))

                if len(row) < len(header):
                    row += [""] * (len(header) - len(row))

                asv = row[asv_col_index].strip()
                if not asv:
                    continue

                rec: Dict[str, str] = {}
                j = 0
                for i in range(len(header)):
                    if i == asv_col_index:
                        continue
                    rec[out_cols[j]] = row[i].strip() if i < len(row) else ""
                    j += 1
                data[asv] = rec

    return out_cols, data


def main():
    ap = argparse.ArgumentParser(description="Merge ASV classification files into one TSV with a 2-line header.")
    ap.add_argument("--input", action="append", required=True,
                    help="SOURCE:PATH:ASV_COLUMN (1-based). Repeat for each file.")
    ap.add_argument("-o", "--output", default="merged_classifications.tsv", help="Output TSV file")
    args = ap.parse_args()

    sources: List[Tuple[str, List[str], Dict[str, Dict[str, str]]]] = []
    all_asvs = set()

    for inp in args.input:
        source, path, asv_idx = parse_input(inp)
        cols, dat = read_table(path, asv_idx)
        sources.append((source, cols, dat))
        all_asvs.update(dat.keys())

    # ordered (source, col) pairs
    column_keys: List[Tuple[str, str]] = []
    for source, cols, _ in sources:
        for c in cols:
            column_keys.append((source, c))

    # quick lookup of per-source data dict
    data_by_source = {source: dat for source, _, dat in sources}

    with open(args.output, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")

        # header row 1: source labels
        h1 = ["ASV_ID"]
        for source, cols, _ in sources:
            h1.extend([source] * len(cols))
        w.writerow(h1)

        # header row 2: original column names
        h2 = ["ASV_ID"]
        for _, cols, _ in sources:
            h2.extend(cols)
        w.writerow(h2)

        # rows
        for asv in sorted(all_asvs, key=asv_sort_key):
            row = [asv]
            for source, col in column_keys:
                row.append(data_by_source[source].get(asv, {}).get(col, ""))
            w.writerow(row)

    print(f"Wrote merged file: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
