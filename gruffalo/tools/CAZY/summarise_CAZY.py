#   #!/usr/bin/python

'''
This Script is intended to summarise the results of CAZY, giving information on
which proteins contain a signal peptide and giving breakdowns of the numbers of
CAZY genes by family and by contig.
'''

import sys
import argparse
import re
from sets import Set
from collections import defaultdict
from collections import Counter
from operator import itemgetter


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

ap = argparse.ArgumentParser()
ap.add_argument('--cazy',required=True,type=str,help='input CAZY file')
ap.add_argument('--inp_secreted',required=True,type=str,help='List of secreted genes')
ap.add_argument('--inp_gff',required=True,type=str,help='gene annotations to provide a breakdown by contig')
ap.add_argument('--summarise_family',required=False,default=False, action='store_true',help='Summarise CAZY gene models by CAZY functional family')
ap.add_argument('--by_contig',required=False,default=False, action='store_true',help='Summarise CAZY gene models by contig rather than accross the whole genome')
ap.add_argument('--trim_gene_id',required=False,default='0', type=int,help='trim the final xbp from each gene id, so that numbers reflect gene numbers rather than proteins')
ap.add_argument('--kubicek_2014',required=False,default=False, action='store_true',help='Summarise CAZY families by functions described in Kubicek et al. 2014')

# ap.add_argument('--out_secreted',required=True,type=str,help='headers of secreted CAZY genes')
# ap.add_argument('--out_by_contig',required=True,type=str,help='breakdown of secreted CAZY gene families by contig')

conf = ap.parse_args() #sys.argv

with open(conf.cazy) as f:
    inp_lines = f.readlines()

with open(conf.inp_secreted) as f:
    secreted_lines = f.readlines()

with open(conf.inp_gff) as f:
    gff_lines = f.readlines()

summarise_family = conf.summarise_family
by_contig = conf.by_contig
trim_id = conf.trim_gene_id
kubicek_2014 = conf.kubicek_2014

# kubicek_2014_dict =
kubicek_2014_dict = defaultdict(str, {
    "GH6.hmm" : "Cellulases",
    "GH7.hmm" : "Cellulases",
    "GH5.hmm" : "Cellulases",
    "GH12.hmm" : "Cellulases",
    "GH45.hmm" : "Cellulases",
    "GH1.hmm" : "B-Glycosidases",
    "GH3.hmm" : "B-Glycosidases",
    "GH61.hmm" : "Accessory enzymes",
    "CBM.hmm" : "Accessory enzymes",
    "GH10.hmm" : "Xylanases",
    "GH11.hmm" : "Xylanases",
    "GH30.hmm" : "Xylanases",
    "GH74.hmm" : "Xyloglucanases",
    "GH27.hmm" : "A-Galactosidases",
    "GH36.hmm" : "A-Galactosidases",
    "GH26.hmm" : "B-Mannase",
    "GH43.hmm" : "A-Arabinosidases",
    "GH51.hmm" : "A-Arabinosidases",
    "GH54.hmm" : "A-Arabinosidases",
    "GH62.hmm" : "A-Arabinosidases",
    "GH35.hmm" : "B-Galactosidases",
    "GH67.hmm" : "B-Glucuronidases",
    "GH115.hmm" : "B-Glucuronidases",
    "GH28.hmm" : "Polygalacturonase",
    "GH78.hmm" : "Polygalacturonase",
    "PL1.hmm" : "Polygalacturonate lyases",
    "PL3.hmm" : "Polygalacturonate lyases",
    "PL4.hmm" : "Polygalacturonate lyases",
    "PL9.hmm" : "Polygalacturonate lyases",
    "PL11.hmm" : "Polygalacturonate lyases"
    }
)


#-----------------------------------------------------
# Step 2
# Create a set of secreted genes
#-----------------------------------------------------

secreted_set = Set([])
for line in secreted_lines:
    line = line.rstrip()
    split_line = line.split()
    gene_id = split_line[0].replace(">", "")
    if trim_id > 0:
        gene_id = gene_id[:-int(trim_id)]
    secreted_set.add(gene_id)

# print secreted_set


#-----------------------------------------------------
# Step 3
# Process CAZY file
#-----------------------------------------------------

# Hmm2function_dict = { }

seen_set = Set([])
family_set = Set([])
gene_dict = defaultdict(list)
family_dict = defaultdict(list)
for line in inp_lines:
    line = line.rstrip()
    split_line = line.split()
    family = split_line[0]
    if summarise_family == True:
        if kubicek_2014 == True:
            # print family
            if kubicek_2014_dict[family]:
                family = kubicek_2014_dict[family]
            else:
                family = "other"
        else:
            family = family[0:2]
            family = family.replace("CB", "CBM")
    gene_id = split_line[2]
    if trim_id > 0:
        gene_id = gene_id[:-int(trim_id)]
    # print gene_id
    if gene_id in secreted_set and gene_id not in seen_set:
        seen_set.add(gene_id)
        family_set.add(family)
        gene_dict[gene_id].append(family)
        family_dict[family].append(gene_id)

for family in family_dict.keys():
    family_sz = len(family_dict[family])
    if by_contig == False:
        print " - ".join([str(family), str(family_sz)])

if by_contig == False:
    quit()

#-----------------------------------------------------
# Step 4
# By contig
#-----------------------------------------------------

gene_count_dict = defaultdict(int)
contig_dict = defaultdict(list)

for line in gff_lines:
    line = line.rstrip()
    split_line = line.split()
    if split_line[2] == "mRNA":
        col9 = split_line[8]
        col9_split = col9.split(";")
        gene_id = col9_split[0]
        gene_id = gene_id.replace("ID=", "")
        if trim_id > 0:
            gene_id = gene_id[:-int(trim_id)]
        contig = split_line[0]
        gene_count_dict[contig] += 1
        if gene_id in seen_set:
            # family = family_dict[gene_id]
            family = ";".join(gene_dict[gene_id])
            # print family
            # print gene_id
            contig_dict[contig].append(family)

keys = contig_dict.keys()
for contig in sorted(keys):
    family_list = contig_dict[contig]
    num_CAZY_genes = len(family_list)
    gene_count = gene_count_dict[contig]
    print " - ".join([str(contig), str(num_CAZY_genes), str(gene_count)])
    # print family_list
    genes_by_family = Counter(family_list)
    # print genes_by_family
    for family in family_set:
        print "\t".join([family, str(genes_by_family[family])])
