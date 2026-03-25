#!/usr/bin/env python3
from pathlib import Path
from typing import List, Tuple, Iterator
import argparse


def fasta_iter(path: Path) -> Iterator[Tuple[str, str]]:
    header = None
    seq_chunks: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append("".join(line.split()))
        if header is not None:
            yield header, "".join(seq_chunks)


def main():
    ap = argparse.ArgumentParser(
        description="Exact substring match: short sequences inside long sequences."
    )
    ap.add_argument("--short", required=True, type=Path,
                    help="FASTA of short sequences (queries).")
    ap.add_argument("--long", required=True, type=Path,
                    help="FASTA of long sequences (targets).")
    ap.add_argument("-o", "--out", default="short_in_long_matches.tsv", type=Path)
    args = ap.parse_args()

    # Load short sequences into memory
    short_seqs = []
    for h, s in fasta_iter(args.short):
        short_seqs.append((h, s.upper()))

    # Write output
    with args.out.open("w") as out:
        out.write("short_id\tshort_len\tlong_id\tmatch_start_0based\tmatch_end_0based\n")

        for long_h, long_seq in fasta_iter(args.long):
            text = long_seq.upper()
            for short_id, short_seq in short_seqs:
                pos = text.find(short_seq)
                if pos != -1:
                    out.write(
                        f"{short_id}\t{len(short_seq)}\t{long_h}\t{pos}\t{pos + len(short_seq) - 1}\n"
                    )

    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
