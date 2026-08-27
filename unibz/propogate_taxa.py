#!/usr/bin/env python3

import pandas as pd
import argparse

# -----------------------------
# Argument parsing
# -----------------------------
parser = argparse.ArgumentParser(
    description="Replace 'unidentified' taxonomy entries with hierarchical 'unclassified_' labels"
)

parser.add_argument(
    "-i", "--input",
    required=True,
    help="Input TSV file"
)

parser.add_argument(
    "-o", "--output",
    default="fixed_taxonomy.tsv",
    help="Output TSV file (default: fixed_taxonomy.tsv)"
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
# Taxonomy columns (edit if needed)
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

# Ensure all expected columns exist
missing = [c for c in tax_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing taxonomy columns: {missing}")

# -----------------------------
# Main function
# -----------------------------
def fix_row(row):
    last_known = None

    for col in tax_cols:
        val = str(row[col]).strip()

        if val.lower() in ["unidentified", "unclassified", "Incertae_Sedis",  "", "nan", "none"]:
            if last_known is not None:
                row[col] = f"unclassified_{last_known}"
            else:
                row[col] = "unclassified_root"
        else:
            last_known = val

    return row

# Apply transformation
df = df.apply(fix_row, axis=1)

# -----------------------------
# Save output
# -----------------------------
df.to_csv(args.output, sep="\t", index=False)

print(f"Saved cleaned taxonomy to: {args.output}")
