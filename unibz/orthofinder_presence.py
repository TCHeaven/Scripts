#!/usr/bin/env python3

import pandas as pd
import argparse


parser = argparse.ArgumentParser()

parser.add_argument(
    "--orthogroups",
    required=True,
    help="Orthogroups.tsv from OrthoFinder"
)

parser.add_argument(
    "--out",
    default="orthogroup_presence.tsv"
)

args = parser.parse_args()


df = pd.read_csv(
    args.orthogroups,
    sep="\t"
)


# first column is OG name
og = df.iloc[:,0]

presence = pd.DataFrame()

presence["Orthogroup"] = og


for col in df.columns[1:]:
    presence[col] = df[col].notna() & (df[col] != "")


presence.to_csv(
    args.out,
    sep="\t",
    index=False
)


print(presence.sum())
