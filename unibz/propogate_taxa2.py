#!/usr/bin/env python3

import pandas as pd
import argparse

# -----------------------------
# Argument parsing
# -----------------------------
parser = argparse.ArgumentParser(
    description="Replace 'unidentified' and 'Incertae_sedis' taxonomy entries with hierarchical 'unclassified_' labels"
)

parser.add_argument(
    "-i", "--input",
    required=True,
    help="Input TSV file"
)

parser.add_argument(
    "-o", "--output",
    default="fixed_taxonomy.tsv",
    help="Output TSV file"
)

parser.add_argument(
    "--sep",
    default="\t",
    help="Column separator (default: tab)"
)

args = parser.parse_args()

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(args.input, sep=args.sep)

# -----------------------------
# Taxonomy columns
# -----------------------------
tax_cols = [
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
]

missing = [c for c in tax_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing taxonomy columns: {missing}")

# -----------------------------
# Helper: identify "bad" taxonomy labels
# -----------------------------
def is_unclassified(val):
    val = str(val).strip().lower()
    
    return (
        val in ["", "nan", "none", "unidentified", "unclassified"]
        or "incertae" in val
    )

# -----------------------------
# Main function
# -----------------------------
def fix_row(row):
    last_known = None

    for col in tax_cols:
        val = row[col]

        if is_unclassified(val):
            if last_known is not None:
                row[col] = f"unclassified_{last_known}"
            else:
                row[col] = "unclassified_root"
        else:
            last_known = val

    return row

# Apply
df = df.apply(fix_row, axis=1)

# -----------------------------
# Save
# -----------------------------
df.to_csv(args.output, sep="\t", index=False)

print(f"Saved cleaned taxonomy to: {args.output}")
