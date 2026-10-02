#!/usr/bin/env python3
import pandas as pd
import argparse
import os
from functools import lru_cache

parser = argparse.ArgumentParser(description="MEGAN-style LCA with rank filtering + full taxonomy path")
parser.add_argument("--input", nargs="+", required=True)
parser.add_argument("--nodes", required=True)
parser.add_argument("--names", required=True)
parser.add_argument("--min-rank", default="genus")
parser.add_argument("--outdir", default="lca_results")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# -----------------------------
# Load taxonomy
# -----------------------------
RANK_ALIAS = {"domain": "superkingdom"}  # newer NCBI dumps use "domain"; remove to keep superkingdom=None

def load_nodes(nodes_file):
    parent, rank_map = {}, {}
    with open(nodes_file) as f:
        for line in f:
            parts = line.split("\t|\t")
            taxid = int(parts[0])
            parent[taxid] = int(parts[1])
            r = parts[2].strip()
            rank_map[taxid] = RANK_ALIAS.get(r, r)
    return parent, rank_map

def load_names(names_file):
    names = {}
    with open(names_file) as f:
        for line in f:
            parts = line.split("\t|\t")
            if "scientific name" in parts[3]:
                names[int(parts[0])] = parts[1]
    return names

parent_map, rank_map = load_nodes(args.nodes)
names_map = load_names(args.names)

# -----------------------------
# Rank system
# -----------------------------
rank_order = ["superkingdom", "kingdom", "phylum", "class",
              "order", "family", "genus", "species"]
rank_index = {r: i for i, r in enumerate(rank_order)}
min_rank_idx = rank_index.get(args.min_rank, 6)

@lru_cache(maxsize=None)
def keep_taxid(taxid):
    """Keep hit if it sits at or below the min rank."""
    current = taxid
    while current in parent_map:
        r = rank_map.get(current, "no rank")
        if r in rank_index and rank_index[r] >= min_rank_idx:
            return True
        if parent_map[current] == current:
            break
        current = parent_map[current]
    return False

# -----------------------------
# LCA helpers
# -----------------------------
@lru_cache(maxsize=None)
def get_ancestors(taxid):
    path = set()
    while taxid in parent_map and taxid != parent_map[taxid]:
        path.add(taxid)
        taxid = parent_map[taxid]
    return frozenset(path)

@lru_cache(maxsize=None)
def get_depth(taxid):
    depth = 0
    while taxid in parent_map and taxid != parent_map[taxid]:
        taxid = parent_map[taxid]
        depth += 1
    return depth

def lca_taxids(taxids):
    paths = [get_ancestors(t) for t in taxids if t is not None]
    if not paths:
        return None
    common = frozenset.intersection(*paths)
    if not common:
        return None
    return max(common, key=get_depth)

# -----------------------------
# MEGAN bit-score filter
# -----------------------------
def filter_hits(df):
    best = df["bit_score"].max()
    return df[df["bit_score"] >= best * 0.90]

# -----------------------------
# Full taxonomy path
# -----------------------------
def build_full_taxonomy(taxid):
    ranks = {r: None for r in rank_order}   # rank -> taxon name
    current = taxid
    while current in parent_map:
        r = rank_map.get(current, "no rank")
        if r in ranks and ranks[r] is None:
            ranks[r] = names_map.get(current, f"taxid_{current}")
        if parent_map[current] == current:
            break
        current = parent_map[current]

    # Propagate: every empty rank BELOW the first assigned rank gets
    # "unclassified_<name of nearest assigned ancestor>"
    last_name = None
    out = {}
    for r in rank_order:
        if ranks[r] is not None:
            last_name = ranks[r]
            out[r] = ranks[r]
        elif last_name is not None:
            out[r] = f"unclassified_{last_name}"
        else:
            out[r] = None
    return out

# -----------------------------
# Process file
# -----------------------------
def process_file(blast_file):
    print(f"\n[INFO] Processing {blast_file}")
    cols = ["read_id", "subject", "pid", "aln_len", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bit_score",
            "taxid", "taxon"]
    df = pd.read_csv(blast_file, sep="\t", names=cols, low_memory=False)
    print("RAW rows:", len(df))
    df["taxid"] = df["taxid"].astype(str).str.split(";").str[0]
    df = df[df["taxid"].str.isnumeric()].copy()
    df["taxid"] = df["taxid"].astype(int)
    df = df.groupby("read_id", group_keys=False).apply(filter_hits)
    print("After bit filter:", len(df))
    if df.empty:
        print("[WARNING] No hits left after bit-score filtering")
        return
    df = df[df["taxid"].apply(keep_taxid)]
    print("After rank filter:", len(df))
    if df.empty:
        print("[WARNING] No hits left after rank filtering")
        return
    lca = df.groupby("read_id")["taxid"].apply(lambda x: lca_taxids(tuple(x.unique())))
    print("LCA NULL:", lca.isna().sum())
    print("Total reads with LCA:", len(lca))
    cache = {}
    rows = []
    for read_id, taxid in lca.items():
        if pd.isna(taxid):
            continue
        taxid = int(taxid)
        if taxid not in cache:
            cache[taxid] = build_full_taxonomy(taxid)
        rows.append(cache[taxid])
    print("ROWS GENERATED:", len(rows))
    print("SAMPLE ROW:", rows[0] if rows else "EMPTY")
    if not rows:
        print("[WARNING] No LCA assignments produced")
        return
    abundance_df = pd.DataFrame(rows)
    # dropna=False is essential: None in superkingdom would otherwise drop every row
    abundance = (abundance_df
                 .groupby(rank_order, dropna=False)
                 .size()
                 .reset_index(name="read_count"))
    abundance["relative_abundance"] = abundance["read_count"] / abundance["read_count"].sum()
    base = os.path.splitext(os.path.basename(blast_file))[0]
    out_file = os.path.join(args.outdir, f"{base}.min-{args.min_rank}.lca.tsv")
    abundance.to_csv(out_file, sep="\t", index=False)
    print(f"[DONE] {out_file}")

for f in args.input:
    process_file(f)
