#!/usr/bin/python


'''
This tool converts gff files into a tab-based format that foramat_embly.py can
use, along with a fasta file to produce an EMBL format file of the genome and
its annotations. THis EMBL file can then be provided to antismash for prediction
of antimsash clusters. This program requires a gff file as input.
'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys,argparse
import re
# from natsort import natsorted, ns
from Bio.SeqFeature import SeqFeature, FeatureLocation
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument('--fasta',required=True,type=str,help='.fasta file of the assembly')
ap.add_argument('--gff',required=True,type=str,help='.gff file of gene annotations')
ap.add_argument('--out',required=True,type=str,help='output directory for files')

conf = ap.parse_args()

outdir = conf.out

#-----------------------------------------------------
# Step 1
# Prepare an individual fasta file for each assembled contig
#-----------------------------------------------------

with open(conf.fasta) as f:
    fasta_lines = f.readlines()
f.close()

First = True
for line in fasta_lines:
    line = line.rstrip()
    # print line
    if line.startswith(">"):
        header = line.replace(">", "")
        if First == True:
            First = False
        else:
            out_file.close()
        out_file = open(outdir + "/" + header + ".fasta","w")
        out_file.write(line + "\n")
    else:
        out_file.write(line + "\n")
out_file.close()

#-----------------------------------------------------
# Step 2
# Parse gff annotations into an annotation format
# with locations of associated fasta files included in
# the first column
#-----------------------------------------------------


with open(conf.gff) as f:
    gff_lines = f.readlines()
f.close()

contig_dict = defaultdict(list)

outlines = []


for line in gff_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    # print split_line
    if split_line[2] and "gene" in split_line[2]:
        # print split_line
        contig = split_line[0] + ".fasta"
        locus_tag = split_line[8].replace("ID=", "").replace(";", "")
        start = int(split_line[3])
        end = int(split_line[4])
        annotation = "hypothetical protein"
        contig_dict[contig].append([contig, locus_tag, start, end, annotation])
        # outlines.append("\t".join([contig, locus_tag, start, end, annotation]))


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


out_file = open(outdir + "/annotationtable.txt","w")
out_file.write("\t".join(["contig", "locus tag", "start", "end", "annotation"]) + "\n")

keys = contig_dict.keys()
for key in natural_sort(keys):
    for outline in sorted(contig_dict[key], key=itemgetter(2)):
        # print key + "\t" + "\t".join(str(x) for x in outline)
        out_file.write("\t".join(str(x) for x in outline) + "\n")
out_file.close()



# print "\n".join(outlines)
