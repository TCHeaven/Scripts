import sys
from collections import defaultdict
import random
import numpy as np
from cyvcf2 import VCF
import pandas as pd
import openpyxl

random.seed(1234)

arg1 = sys.argv[1]
arg2 = sys.argv[2]
arg3 = sys.argv[3]
arg4 = sys.argv[4]
arg5 = sys.argv[5]
arg6 = sys.argv[6]
arg7 = sys.argv[7]
print(arg1)
print(arg2)
print(arg3)
print(arg4)
print(arg5)
print(arg6)
print(arg7)

vcf_file = arg1

samples = arg2.split(",")

chroms = arg3.split(",")

scaffold_lengths = {}
scaffold_sizes_list = arg4.split(",")  # Assuming scaffold lengths are provided as comma-separated values

for scaffold_size in scaffold_sizes_list:
    scaffold, length = scaffold_size.split(":")
    scaffold_lengths[scaffold.strip()] = int(length.strip())


iterations = arg5.split(",")

replicates = int(arg6)

OutFile = arg7

allele_dict = {
    0: 'HOM_REF',
    1: 'HET',
    3: 'UNKNOWN',
    2: 'HOM_VAR',
}

vcf = VCF(vcf_file, gts012=True)

target_index = [vcf.samples.index(x) for x in samples]  
target_index2samples = {vcf.samples.index(x): x for x in samples}  



def get_all_covered_genome_positions(vcf_object, chromsomes, scaffold_lengths):
    """

    """

    univ = set()
    box_addresses = defaultdict(set)
    # load all the ticks
    for chrom in chroms:  # EVERY Chromosome
        scaffold_length = scaffold_lengths[chrom]
        for v in vcf_object(chrom + ':' + str(0) + '-' + str(scaffold_lengths[chrom])):  # every positon
            target_gts_composition = {samples[i]: allele_dict[v.gt_types[i]] for i in target_index}
            for this_sample, letter in target_gts_composition.items():
                if letter not in ['HOM_VAR']:
                    continue
                address = str(v.CHROM) + '@' + str(v.POS)  # every sample
                box_addresses[this_sample].add(address)
                univ.add(address)

    universal_ticks_len = len(univ)

    return box_addresses, universal_ticks_len


if __name__ == "__main__":

    box, universal_ticks = get_all_covered_genome_positions(vcf, chroms, scaffold_lengths)
    print(f'Number of unique positions covered in the ref genome = {universal_ticks}')  # 100 percent

    with open(f'{OutFile}-samples.csv', 'w') as outf:
        outf.write(','.join([str(x) for x in ['sample', 'no_of_unique_polymorphisms']]) + '\n')
        for _, sample in target_index2samples.items():
            if sample in box:
                outf.write(','.join([str(x) for x in [sample, len(box[sample])]]) + '\n')

    with open(f'{OutFile}-hits.tab', 'w') as outf, open(f'{OutFile}-curve.csv',
                                                                         'w') as outC:
        outf.write('\t'.join([str(x) for x in ['sample_size', 'replicate', 'uniq_polymorphisms', 'boot_samples']]) + '\n')
        outC.write(','.join([str(x) for x in ['sample_size', 'rep_polymorphisms_median', 'rep_polymorphisms_median_fraction']]) + '\n')

        for sample_size in iterations:
            rep_polymorphisms_this_size = []
            for replicate in range(1, replicates):  # 100 replicates dont needs reps for the full size i.e. 209
                sampled_target_index = random.sample(target_index, int(sample_size))  # random sampling without replacement
                boot_samples = [target_index2samples[i] for i in sampled_target_index]
                assert 'ligustri' not in boot_samples
                # 6 chrom
                seen = set()  # all possible covered positions in a genome # 17.4M ( 100%)
                for sample in boot_samples:
                    seen.update(box[sample])
                uniq_polymorphisms = len(seen)
                print(f' {replicate} {sample_size}  {uniq_polymorphisms} {boot_samples}')  #
                rep_polymorphisms_this_size.append(uniq_polymorphisms)
                outf.write('\t'.join([str(x) for x in [sample_size, replicate, uniq_polymorphisms] + boot_samples]) + '\n')
            # -----------------
            rep_polymorphisms_avg = np.median(rep_polymorphisms_this_size)

            print(f'rep_polymorphisms_median = {rep_polymorphisms_avg}')
            outC.write(','.join([str(x) for x in [sample_size, rep_polymorphisms_avg, rep_polymorphisms_avg / universal_ticks]]) + '\n')

    print('Done')

 
