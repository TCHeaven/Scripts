#!/usr/bin/env python3

import argparse
import pandas as pd
import sys


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--snps", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--max-missing", type=float, default=1.0)
    return parser.parse_args()


def main():
    args = parse_args()
    print("Reading:", args.snps)
    # CRITICAL FIX: use tab, not whitespace
    df = pd.read_csv(
        args.snps,
        sep="\t",
        dtype=str
    )
    print("\nColumns detected:", len(df.columns))
    print(df.columns.tolist())
    if df.shape[1] < 3:
        sys.exit("ERROR: SNP file parsing failed.")
    reference = df.columns[2]
    samples = list(df.columns[2:])
    print("\nReference:", reference)
    print("\nSamples:")
    for s in samples:
        print(" ", s)
    print("\nSNPs before filtering:", len(df))
    genotype = pd.DataFrame(index=df.index)
    for sample in samples:
        vals = []
        for ref, allele in zip(df[reference], df[sample]):
            # Missing values now preserved correctly
            if pd.isna(allele):
                vals.append(pd.NA)
            elif allele == ref:
                vals.append(0)
            else:
                vals.append(1)
        genotype[sample] = vals
    # Missing filter
    missing_fraction = genotype.isna().mean(axis=1)
    genotype = genotype.loc[missing_fraction <= args.max_missing]
    print("After missing filter:", len(genotype))
    # Remove invariant
    genotype = genotype.loc[
        genotype.nunique(axis=1, dropna=True) > 1
    ]
    print("After invariant removal:", len(genotype))
    # SNP IDs
    snp_ids = (
        df.loc[genotype.index, "refSequence"]
        + "_"
        + df.loc[genotype.index, "refPosition"]
    )
    genotype.index = snp_ids
    # Force integer + NA
    genotype = genotype.astype("Int64")
    genotype.to_csv(
        args.out,
        sep="\t",
        index=True,
        index_label="SNP",
        na_rep="NA"
    )
    print("\nWritten:", args.out)
    print("Dimensions:", genotype.shape)


if __name__ == "__main__":
    main()
