import pandas as pd
import argparse
import os

# -----------------------------
# Argument parser
# -----------------------------
parser = argparse.ArgumentParser(description="MEGAN-style LCA with rank filtering")
parser.add_argument("--input", nargs="+", required=True, help="BLAST files")
parser.add_argument("--nodes", required=True, help="nodes.dmp")
parser.add_argument("--names", required=True, help="names.dmp")
parser.add_argument("--min-rank", default="genus", help="Minimum rank (e.g. genus, family)")
parser.add_argument("--outdir", default="lca_results")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# -----------------------------
# Load taxonomy
# -----------------------------
def load_nodes(nodes_file):
    parent = {}
    rank = {}
    with open(nodes_file) as f:
        for line in f:
            parts = line.split("\t|\t")
            taxid = int(parts[0])
            parent_taxid = int(parts[1])
            tax_rank = parts[2].strip()
            parent[taxid] = parent_taxid
            rank[taxid] = tax_rank
    return parent, rank

def load_names(names_file):
    names = {}
    with open(names_file) as f:
        for line in f:
            parts = line.split("\t|\t")
            taxid = int(parts[0])
            name = parts[1]
            name_class = parts[3]
            if "scientific name" in name_class:
                names[taxid] = name
    return names

parent_map, rank_map = load_nodes(args.nodes)
names_map = load_names(args.names)

# -----------------------------
# Rank order hierarchy
# -----------------------------
rank_order = [
    "superkingdom", "kingdom", "phylum", "class",
    "order", "family", "genus", "species"
]

rank_index = {r: i for i, r in enumerate(rank_order)}
min_rank_idx = rank_index.get(args.min_rank, 6)  # default genus

# -----------------------------
# Get ancestor at minimum rank
# -----------------------------
def lift_to_rank(taxid):
    while taxid in parent_map:
        r = rank_map.get(taxid, "no rank")
        if r in rank_index and rank_index[r] >= min_rank_idx:
            return taxid
        parent = parent_map[taxid]
        if parent == taxid:
            break
        taxid = parent
    return None  # discard if never reaches rank

# -----------------------------
# Get all ancestors
# -----------------------------
def get_ancestors(taxid):
    path = set()
    while taxid in parent_map and taxid != parent_map[taxid]:
        path.add(taxid)
        taxid = parent_map[taxid]
    return path

# -----------------------------
# LCA
# -----------------------------
def lca_taxids(taxids):
    paths = [get_ancestors(t) for t in taxids if t is not None]
    if not paths:
        return None
    common = set.intersection(*paths)
    return max(common) if common else None

# -----------------------------
# Bit-score filter
# -----------------------------
def filter_hits(group):
    best = group["bit_score"].max()
    return group[group["bit_score"] >= best * 0.9]

# -----------------------------
# Process one file
# -----------------------------
def process_file(blast_file):
    print(f"[INFO] Processing {blast_file}")
    cols = [
        "read_id", "subject", "pid", "aln_len",
        "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send",
        "evalue", "bit_score",
        "taxid", "taxon"
    ]
    df = pd.read_csv(blast_file, sep="\t", names=cols)
    # clean taxid
    df["taxid"] = df["taxid"].astype(str).str.split(";").str[0]
    df = df[df["taxid"].str.isnumeric()]
    df["taxid"] = df["taxid"].astype(int)
    # filter by bit score
    df = df.groupby("read_id", group_keys=False).apply(filter_hits)
    # lift to minimum rank
    df["taxid"] = df["taxid"].apply(lift_to_rank)
    # drop bad ones
    df = df.dropna(subset=["taxid"])
    # LCA per read
    lca = df.groupby("read_id")["taxid"].apply(lambda x: lca_taxids(x.unique()))
    # map to names
    lca_names = lca.map(lambda x: names_map.get(x, f"taxid_{x}") if x else "unclassified")
    # abundance
    abundance = lca_names.value_counts().reset_index()
    abundance.columns = ["taxon", "read_count"]
    abundance["relative_abundance"] = abundance["read_count"] / abundance["read_count"].sum()
    # save
    base = os.path.basename(blast_file).replace(".out", "")
    out_file = os.path.join(
        args.outdir,
        f"{base}.min-{args.min_rank}.lca.tsv"
    )
    abundance.to_csv(out_file, sep="\t", index=False)
    print(f"[DONE] {out_file}")

# -----------------------------
# Run all files
# -----------------------------
for f in args.input:
    process_file(f)
