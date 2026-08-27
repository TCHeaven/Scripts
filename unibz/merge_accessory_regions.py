#!/usr/bin/env python3

import pandas as pd
import argparse


parser = argparse.ArgumentParser(
    description="Merge accessory regions using presence/absence matrix"
)

parser.add_argument("--haplotypes", required=True)
parser.add_argument("--accessory", required=True)
parser.add_argument("--output", required=True)

args = parser.parse_args()


# accessory node indices (as integers)
with open(args.accessory) as f:
    accessory = set(int(x.strip()) for x in f)


df = pd.read_csv(args.haplotypes, sep="\t")


node_cols = [c for c in df.columns if c.startswith("node.")]

# map column -> node index
col_to_node = {
    col: int(col.replace("node.", ""))
    for col in node_cols
}


regions = []


for _, row in df.iterrows():
    path = row["path.name"]
    current = []
    for col in node_cols:
        node_id = col_to_node[col]
        val = row[col]
        if val == 1 and node_id in accessory:
            current.append(node_id)
        else:
            if current:
                regions.append([
                    path,
                    current[0],
                    current[-1],
                    len(current)
                ])
                current = []
    if current:
        regions.append([
            path,
            current[0],
            current[-1],
            len(current)
        ])


out = pd.DataFrame(
    regions,
    columns=["path", "start_node", "end_node", "node_count"]
)

out.to_csv(args.output, sep="\t", index=False)

print("Merged regions:", len(out))
