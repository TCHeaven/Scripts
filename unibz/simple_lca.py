#!/usr/bin/env python3

import os
import argparse
import pandas as pd


# -----------------------------
# 1. Load taxonomy
# -----------------------------
def load_nodes(nodes_file):
    parent = {}
    with open(nodes_file, "r") as f:
        for line in f:
            parts = line.split("\t|\t")
            taxid = int(parts[0])
            parent_taxid = int(parts[1])
            parent[taxid] = parent_taxid
    return parent


def load_names(names_file):
    names = {}
    with open(names_file, "r") as f:
        for line in f:
            parts = line.split("\t|\t")
            taxid = int(parts[0])
            name = parts[1]
            name_class = parts[3]
            if "scientific name" in name_class:
                names[taxid] = name
    return names


# -----------------------------
# 2. Taxonomy helpers
# -----------------------------
def get_ancestors(taxid, parent_map):
    path = set()
    while taxid in parent_map and taxid != parent_map[taxid]:
        path.add(taxid)
        taxid = parent_map[taxid]
    return path


def lca_taxids(taxids, parent_map):
    paths = [get_ancestors(t, parent_map) for t in taxids]
    if not paths:
        return None
    common = set.intersection(*paths)
    return max(common) if common else 1


# -----------------------------
# 3. MEGAN-style filtering
# -----------------------------
def filter_hits(df):
    best = df["bit_score"].max()
    threshold = best * 0.90
    return df[df["bit_score"] >= threshold]


# -----------------------------
# 4. Process one file
# -----------------------------
def process_file(blast_file, parent_map, names_map, output_dir):
    print(f"[INFO] Processing: {blast_file}")
    cols = [
        "read_id", "subject", "pid", "aln_len",
        "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send",
        "evalue", "bit_score",
        "taxid", "taxon"
    ]
    df = pd.read_csv(blast_file, sep="\t", names=cols)
    # fix multi-taxid issue
    df["taxid"] = (
        df["taxid"]
        .astype(str)
        .str.split(";")
        .str[0]
        .astype(int)
    )
    # filter hits (MEGAN-style)
    df_filt = df.groupby("read_id", group_keys=False).apply(filter_hits)
    # LCA per read
    lca_series = df_filt.groupby("read_id")["taxid"].apply(
        lambda x: lca_taxids(x.unique(), parent_map)
    )
    # map to names
    lca_names = lca_series.map(lambda x: names_map.get(x, f"taxid_{x}"))
    # abundance table
    abundance = lca_names.value_counts().reset_index()
    abundance.columns = ["taxon", "read_count"]
    abundance["relative_abundance"] = (
        abundance["read_count"] / abundance["read_count"].sum()
    )
    # output file name
    base = os.path.basename(blast_file).replace(".out", "")
    out_file = os.path.join(output_dir, f"{base}_lca_abundance.tsv")
    abundance.to_csv(out_file, sep="\t", index=False)
    print(f"[DONE] Saved: {out_file}")
    print(abundance.head(10))
    return abundance


# -----------------------------
# 5. Main
# -----------------------------
def main():
    parser = argparse.ArgumentParser(description="MEGAN-style LCA pipeline (batch mode)")
    parser.add_argument("--input", nargs="+", required=True, help="BLAST output files")
    parser.add_argument("--nodes", required=True, help="nodes.dmp file")
    parser.add_argument("--names", required=True, help="names.dmp file")
    parser.add_argument("--outdir", default="lca_results", help="output directory")
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    print("[INFO] Loading taxonomy...")
    parent_map = load_nodes(args.nodes)
    names_map = load_names(args.names)
    for f in args.input:
        process_file(f, parent_map, names_map, args.outdir)


if __name__ == "__main__":
    main()
