#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import Iterator, Tuple, List, Dict


def fasta_iter(path: Path) -> Iterator[Tuple[str, str]]:
    """
    Stream FASTA records from `path`.
    Returns (header_without_>, sequence_with_no_whitespace).
    Handles multi-line sequences.
    """
    header = None
    seq_chunks: List[str] = []

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                # remove all whitespace in sequence lines
                seq_chunks.append("".join(line.split()))
        if header is not None:
            yield header, "".join(seq_chunks)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Deduplicate identical reads in a FASTA file and emit unique FASTA + mapping TSV."
    )
    ap.add_argument("input_fasta", type=Path, help="Input FASTA (can be multi-line sequences).")
    ap.add_argument("-o", "--out-fasta", type=Path, default=None,
                    help="Output FASTA with unique sequences. Default: <input>.unique.fasta")
    ap.add_argument("-m", "--out-map", type=Path, default=None,
                    help="Output TSV mapping. Default: <input>.unique.map.tsv")
    ap.add_argument("--prefix", default="uniq", help="Prefix for unique IDs (default: uniq).")
    ap.add_argument("--id-width", type=int, default=7, help="Zero-padding width for IDs (default: 7).")
    args = ap.parse_args()

    in_path: Path = args.input_fasta
    if args.out_fasta is None:
        out_fa = in_path.with_suffix(in_path.suffix + ".unique.fasta")
    else:
        out_fa = args.out_fasta

    if args.out_map is None:
        out_map = in_path.with_suffix(in_path.suffix + ".unique.map.tsv")
    else:
        out_map = args.out_map

    # seq -> unique index (1-based)
    seq2idx: Dict[str, int] = {}
    rep_header: Dict[int, str] = {}
    members: Dict[int, List[str]] = {}
    counts: Dict[int, int] = {}

    n_unique = 0
    n_total = 0

    for hdr, seq in fasta_iter(in_path):
        n_total += 1
        if seq in seq2idx:
            idx = seq2idx[seq]
            counts[idx] += 1
            members[idx].append(hdr)
        else:
            n_unique += 1
            idx = n_unique
            seq2idx[seq] = idx
            rep_header[idx] = hdr
            members[idx] = [hdr]
            counts[idx] = 1

    # Write unique FASTA
    with out_fa.open("w", encoding="utf-8") as fa_out:
        # We can reconstruct sequences by reversing seq2idx (store a list)
        # For memory simplicity: build list indexed by idx
        uniq_seqs: List[str] = [""] * (n_unique + 1)
        for seq, idx in seq2idx.items():
            uniq_seqs[idx] = seq

        for idx in range(1, n_unique + 1):
            uid = f"{args.prefix}{idx:0{args.id_width}d}"
            fa_out.write(f">{uid} count={counts[idx]} rep={rep_header[idx]}\n")
            fa_out.write(uniq_seqs[idx] + "\n")

    # Write mapping TSV (one member per row)
    with out_map.open("w", encoding="utf-8") as map_out:
        map_out.write("unique_id\trepresentative_header\tmember_header\n")
        for idx in range(1, n_unique + 1):
            uid = f"{args.prefix}{idx:0{args.id_width}d}"
            rep = rep_header[idx]
            for mh in members[idx]:
                map_out.write(f"{uid}\t{rep}\t{mh}\n")

    print(f"Input reads: {n_total}")
    print(f"Unique reads: {n_unique}")
    print(f"Wrote: {out_fa}")
    print(f"Wrote: {out_map}")


if __name__ == "__main__":
    main()
