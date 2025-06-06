#!/usr/bin/python


'''
Tool to extract lengths of 5' and 3' te lengths for all genes in the
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

#Provide a contcatenated .gff file with gene and Earlgrey annotations
ap.add_argument('--Gff',required=True,type=str,help='Gene models in Gff format')

conf = ap.parse_args()

with open(conf.Gff) as f:
    gff_lines = f.readlines()

# set a dictionary to hold gene information for each contig
te_dict = defaultdict(list)
contig_set = Set()
contig_list = []


# for each gene extract start, stop, gene name and orientation
# for each te extract start, stop, family, ID and orientation informations
# sort start and stop by size
# store the gene information in a dictionary as a list of lists [start, stop, gene, orientation]
for line in gff_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    tool = split_line[1]
    feature = split_line[2]
    if tool == 'Earl_Grey':
        contig = split_line[0]
        if contig not in contig_set:
            contig_set.add(contig)
            contig_list.append(contig)
        start = split_line[3]
        end = split_line[4]
        strand = split_line[6]
        col_3 = split_line[2]
        Super_ID = col_3.split("/")[0].split()[0]
        ID = col_3.split()[0]
        te_dict[contig].append([start, end, strand, Super_ID, ID])
    elif feature == 'gene':
        contig = split_line[0]
        if contig not in contig_set:
            contig_set.add(contig)
            contig_list.append(contig)
        start = split_line[3]
        end = split_line[4]
        strand = split_line[6]
        col_9 = split_line[8]
        Super_ID = 'Gene'
        ID = col_9.split(";")[0]
        ID = ID.replace("ID=", "")
        te_dict[contig].append([start, end, strand, Super_ID, ID])

# Once all genes have been entered into the dictionary, for each key sort the gene entries
# upon start location



for contig in contig_list:
    # print contig
    # print gene_dict[key]
    te_list = te_dict[contig]
    te_list.sort(key = lambda x: int(x[0]))
    # print gene_list
    # print len(gene_list)
    for i, element in enumerate(te_list):
        start = element[0]
        end = element[1]
        strand = element[2]
        Super_ID = element[3]
        ID = element[4]
        if Super_ID != 'Gene':
            continue
        # print i
        if i > 0:
            upstream_element = te_list[(i-1)]
            if upstream_element[3] == 'Gene':
                upstream_lgth = 'gene'
            else:
                upstream_lgth = int(start) - int(upstream_element[1])
                if 0 >= upstream_lgth:
                    upstream_lgth = 'Nested TE'
        else:
            upstream_lgth = 'end'
            upstream_element = None
        if i + 1 < len(te_list):
            downstream_element = te_list[(i+1)]
            if downstream_element[3] == 'Gene':
                downstream_lgth = 'gene'
            else:
                downstream_lgth = int(downstream_element[0]) - int(end)
                if 0 >= downstream_lgth:
                    downstream_lgth = 'Nested TE'
        else:
            downstream_lgth = 'end'
            downstream_element = None
        #---
        if strand == '+':
            five_prime_lgth = str(upstream_lgth)
            if upstream_element is None:
                five_prime_family = 'end'
                five_prime_ID = 'end'
		five_prime_strand = 'end'
            else:
                five_prime_family = str(upstream_element[3])
                five_prime_ID = str(upstream_element[4])
		five_prime_strand = str(upstream_element[2])
            three_prime_lgth = str(downstream_lgth)
            if downstream_element is None:
                three_prime_family = 'end'
                three_prime_ID = 'end'
		three_prime_strand = 'end'
            else:
                three_prime_family = str(downstream_element[3])
                three_prime_ID = str(downstream_element[4])
		three_prime_strand = str(downstream_element[2])
        elif strand == '-':
            five_prime_lgth = str(downstream_lgth)
            if downstream_element is None:
                five_prime_family = 'end'
                five_prime_ID = 'end'
		five_prime_strand = 'end'
            else:
                five_prime_ID = str(downstream_element[4])
                five_prime_family = str(downstream_element[3])
		five_prime_strand = str(downstream_element[2])
            three_prime_lgth = str(upstream_lgth)
            if upstream_element is None:
                three_prime_family = 'end'
                three_prime_ID = 'end'
		three_prime_strand = 'end'
            else:
                three_prime_family = str(upstream_element[3])
                three_prime_ID = str(upstream_element[4])
		three_prime_strand = str(upstream_element[2])
        print "\t".join([ID, five_prime_lgth, five_prime_family, five_prime_ID, five_prime_strand, three_prime_lgth, three_prime_family, three_prime_ID, three_prime_strand, strand])

    # exit()

# Itterate through the list of genes on each contig. Store the stop coordinate and
# determine the distance to the start of the next gene.
# stote the upstream and downstream distance of each gene, with 5' and 3' distance
# determined by orientation.

