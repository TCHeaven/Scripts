#!/usr/bin/env python3

import sys
import math

def usage():
    return "\nUSAGE: {} fastaFile windowSize stepSize chromosomeName[forCircos] posColor[forCircos] negColor[forCircos]\n\n".format(sys.argv[0])

def get_fasta_names_and_seqs(filename):
    fasta_names = []
    fasta_seqs = []
    with open(filename, 'r') as file:
        seq = ""
        for line in file:
            line = line.strip()
            if not line:
                continue
            elif line.startswith('>'):
                if seq:
                    fasta_seqs.append(seq.replace(' ', ''))
                fasta_names.append(line[1:])
                seq = ""
            else:
                seq += line
        if seq:
            fasta_seqs.append(seq.replace(' ', ''))
    return fasta_names, fasta_seqs

def calculate_gc_skew(seq, window_size, step_size, chromosome_name, pos_color, neg_color):
    length = len(seq)
    left_size = math.ceil((window_size - 1) / 2)
    right_size = math.floor((window_size - 1) / 2)
    results = []
    for i in range(0, length, step_size):
        gcount = 0
        ccount = 0
        left = i - left_size
        right = i + right_size
        if left < 0:
            left = length + left
        if right >= length:
            right = right - length
        if left < right:
            substr = seq[left:right + 1].upper()
            ccount += substr.count('C')
            gcount += substr.count('G')
        else:
            substr1 = seq[0:right + 1].upper()
            substr2 = seq[left:].upper()
            ccount += substr1.count('C') + substr2.count('C')
            gcount += substr1.count('G') + substr2.count('G')
        if gcount + ccount == 0:
            gc_skew = 0
        else:
            gc_skew = (gcount - ccount) / (gcount + ccount)
        color = pos_color if gc_skew > 0 else neg_color
        results.append(f"{chromosome_name}\t{i + 1}\t{i + 1}\t{gc_skew:.4f}\tfill_color={color}")
    return results

def main():
    if len(sys.argv) != 6:
        print(usage())
        sys.exit(1)
    file = sys.argv[1]
    window_size = int(sys.argv[2])
    step_size = int(sys.argv[3])
    pos_color = sys.argv[4]
    neg_color = sys.argv[5]
    fasta_names, fasta_seqs = get_fasta_names_and_seqs(file)
    if not fasta_seqs:
        print("No sequences found in the FASTA file.")
        sys.exit(1)
    for chromosome_name, seq in zip(fasta_names, fasta_seqs):
        results = calculate_gc_skew(seq, window_size, step_size, chromosome_name, pos_color, neg_color)
        for result in results:
            print(result)

if __name__ == "__main__":
    main()
