#   !/usr/bin/env python

'''
sscp_filter.py

This script aids identification of small-secreted cysteine rich proteins.
Predicted secreted proteins are required as input. From these, those that are
shorter than 150 aa, and hhave a cysteine contenet above the input value are
considered to be short and cysteine-rich.

All results are printed in the header of the output fasta file.
'''

import os
import sys
from Bio import SeqIO
import sys,argparse

ap = argparse.ArgumentParser()
ap.add_argument('--inp_fasta',required=True,type=str,help='.fasta file containing query proteins')
ap.add_argument('--max_length',required=True,type=int,help='maximum length of a "small" protein')
ap.add_argument('--threshold',required=True,type=int,help='% threshold above which proteins are considered cysteine rich')
ap.add_argument('--out_fasta',required=True,type=str,help='output results .fasta file of all proteins')

conf = ap.parse_args()
out_records = []
i = 0

print "% cysteine content threshold set to:\t" + str(conf.threshold)
print "maximum length set to:\t" + str(conf.max_length)

for record in SeqIO.parse(conf.inp_fasta, "fasta"):
        seq = record.seq
        seq_len = len(record)
        cys = (( seq.count('C') / float(seq_len) ) * 100)
        record.description += "\t--cysteine%=\t"
        record.description += str("{0:.0f}".format(cys))
        if seq_len <= conf.max_length and cys >= conf.threshold:
                record.description += "\t--SSCP=\tYes\t"
                i += 1
        else:
                record.description += "\t--SSCP=\tNo\t"
        out_records.append(record)

print "No. short-cysteine rich proteins in input fasta:\t" + str(i)
SeqIO.write(out_records, conf.out_fasta, "fasta")



