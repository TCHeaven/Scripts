#!/usr/bin/env python

""" Find PIC and Nei diversity dfrom VCF 

Nei gene diversity definition Kisman, 2003, Plant Pathology

  https://bsppjournals.onlinelibrary.wiley.com/doi/epdf/10.1046/j.1365-3059.2003.00923.x

Nei measure of the average gene diversity per locus(HS)

Polymorphic Information Content(PIC) is first found using 
 https://link.springer.com/article/10.1007/s10528-012-9509-1 
 
As a test we tried the following in PIC_CALC.exe as an input: 
locus1 2  0.05 0.95
locus2 2  0.1  0.9
locus3 2  0.2  0.8
locus4 2  0.3  0.7
locus5 2  0.4  0.6

>> We get
locus1  2  0.0904875000000001
locus2  2  0.1638
locus3  2  0.2688
locus4  2  0.3318
locus5  2  0.3648


"""

__author__ = "Jitender"
__copyright__ = "Copyright 2022, JIC"

import sys
import click
import numpy as np
from collections import defaultdict
from cyvcf2 import VCF

# Nei diversity equation 13  
def get_hs(qs):
     return 1 - qs*qs  - ((1-qs)*(1-qs))    

def pic(p):
    """
    https://link.springer.com/article/10.1007/s10528-012-9509-1
    """ 
    q = 1 - p
    sum_p2 = p*p + q*q
    p_diag = 2.0*p*p*q*q

    return 1 - sum_p2 - p_diag


@click.command(help="provide input VCF file  ")
@click.option('-i', '--infile',    help="Input VCF file")
@click.option('-o', '--outfile',   help="Output tab file name")
def chrom(infile, outfile):
    """
    This script extract PIC scores
    """

    desired_chroms = ['scaffold_1', 'scaffold_2', 'scaffold_3', 'scaffold_4', 'scaffold_5', 'scaffold_6']

    if infile:
       vcf = infile  
       freqD = defaultdict(list)
       for p in VCF(vcf): # or VCF('some.bcf')
           if p.CHROM in desired_chroms: 
              freqD[p.CHROM].append(p.INFO.get('AF'))

       with open(outfile,'w') as outC:
            pipe = ['chrom', 'num_snps', 'MAF', 'Nei', 'PIC']
            outC.write('\t'.join([str(x) for x in pipe]) + '\n')
            for chrom, vect in freqD.items():
                maf = [af if af < 0.5 else 1.0 - af  for af in vect]
                pipe = [chrom, len(vect)] +['{0:.3f}'.format(np.mean(maf))] + [ '{0:.3f}'.format(np.mean([ get_hs(q) for q in vect]))] +  [ '{0:.3f}'.format( np.mean([ pic(q) for q in vect]))]
                outC.write('\t'.join([str(x) for x in pipe]) + '\n')

if __name__ == '__main__':
    exit(chrom())


 