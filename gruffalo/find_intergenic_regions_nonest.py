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
import math
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
            upstream_lgth = int(start) - int(upstream_element[1])
            upstream_orientation = upstream_element[2]
        else:
            upstream_lgth = 99999
            upstream_orientation = None
        if i + 1 < len(gene_list):
            downstream_element = gene_list[(i+1)]
            downstream_lgth = int(downstream_element[0]) - int(end)
            downstream_orientation = downstream_element[2]
        else:
            downstream_lgth = 99999
            downstream_orientation =  None
        #----
        # This section is to deal with a situations where a predicted effector
        # from an ORF is nested within another gene model
        if 0 >= upstream_lgth and i > 0:
            upstream_lgth = None
            upstream_orientation = None
        if 0 >= downstream_lgth and i + 1 < len(gene_list):
            downstream_lgth = None
            downstream_orientation =  None
        #----
        if strand == '+':
                five_prime_lgth = str(upstream_lgth) if upstream_lgth is not None else 'NA'
                three_prime_lgth = str(downstream_lgth) if downstream_lgth is not None else 'NA'
                upstream_strand = upstream_orientation if upstream_orientation is not None else 'NA'
                downstream_strand = downstream_orientation if downstream_orientation is not None else 'NA'
        elif strand == '-':
                five_prime_lgth = str(downstream_lgth) if downstream_lgth is not None else 'NA'
                three_prime_lgth = str(upstream_lgth) if upstream_lgth is not None else 'NA'
                upstream_strand = downstream_orientation if downstream_orientation is not None else 'NA'
                downstream_strand = upstream_orientation if upstream_orientation is not None else 'NA'
        print "\t".join([ID, five_prime_lgth, three_prime_lgth, strand, upstream_strand, downstream_strand ])

#First and last genes in a contig are given 3' or 5' distance of 99,999.
