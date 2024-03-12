#!/usr/bin/python


'''
Tool to extract lengths of 5' and 3' intergenic regions for all genes in the
genome. Genes flanked by a contig break can be chosen to be retained or excluded.
'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()

ap.add_argument('--Gff',required=True,type=str,help='Gene models in Gff format')
ap.add_argument('--keep_break',required=False,action="store_true",help='Set if intergenic lengths should be counted for genes flanking a contig break')


conf = ap.parse_args()

with open(conf.Gff) as f:
    gff_lines = f.readlines()

kb = False
if conf.keep_break:
    kb = True


# set a dictionary to hold gene information for each contig
gene_dict = defaultdict(list)
contig_set = Set()
contig_list = []

# for each gene extract start, stop, gene name and orientation
# sort start and stop by size
# store the gene information in a dictionary as a list of lists [start, stop, gene, orientation]
for line in gff_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    feature = split_line[2]
    if feature == 'gene':
        contig = split_line[0]
        if contig not in contig_set:
            contig_set.add(contig)
            contig_list.append(contig)
        start = split_line[3]
        end = split_line[4]
        strand = split_line[6]
        col_9 = split_line[8]
        ID = col_9.split(";")[0]
        ID = ID.replace("ID=", "")
        gene_dict[contig].append([start, end, strand, ID])

# Once all genes have been entered into the dictionary, for each key sort the gene entries
# upon start location



for contig in contig_list:
    # print contig
    # print gene_dict[key]
    gene_list = gene_dict[contig]
    gene_list.sort(key = lambda x: int(x[0]))
    # print gene_list
    # print len(gene_list)
    for i, element in enumerate(gene_list):
        start = element[0]
        end = element[1]
        strand = element[2]
        ID = element[3]
        # print i
        if i > 0:
            upstream_element = gene_list[(i-1)]
        else:
            continue
        if i + 1 < len(gene_list):
            downstream_element = gene_list[(i+1)]
        else:
            continue
        upstream_lgth = int(start) - int(upstream_element[1])
        downstream_lgth = int(downstream_element[0]) - int(end)
        #----
        # This section is to deal with a situations where a predicted effector
        # from an ORF is nested within another gene model
        if 0 >= upstream_lgth and i > 1:
            upstream_element = gene_list[(i-2)]
            upstream_lgth = int(start) - int(upstream_element[1])
        elif 0 >= upstream_lgth and i <= 1:
            continue
        if 0 >= downstream_lgth and i + 2 < len(gene_list):
            downstream_element = gene_list[(i+2)]
            downstream_lgth = int(downstream_element[0]) - int(end)
        elif 0 >= downstream_lgth and i + 2 >= len(gene_list):
            continue
        #---
        if strand == '+':
            five_prime_lgth = str(upstream_lgth)
            three_prime_lgth = str(downstream_lgth)
        elif strand == '-':
            five_prime_lgth = str(downstream_lgth)
            three_prime_lgth = str(upstream_lgth)
        print "\t".join([ID, five_prime_lgth, three_prime_lgth, strand])

    # exit()

# Itterate through the list of genes on each contig. Store the stop coordinate and
# determine the distance to the start of the next gene.
# stote the upstream and downstream distance of each gene, with 5' and 3' distance
# determined by orientation.
