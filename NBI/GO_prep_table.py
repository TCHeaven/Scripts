#!/usr/bin/python

import sys
import argparse
import re
import numpy
from sets import Set
from Bio import SeqIO
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument('--interpro',required=True,type=str,help='The genome assembly')


conf = ap.parse_args()

with open(conf.interpro) as f:
    interpro_lines = f.readlines()

annot_dict = defaultdict(list)

for line in interpro_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[0]
    m = re.findall("GO:.......", line)
    if m:
        annot_dict[gene_id].extend(m)

for gene_id in annot_dict.keys():
    annotation_set = set(annot_dict[gene_id])
    print gene_id + "\t" + ", ".join(annotation_set)
