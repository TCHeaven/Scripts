import sys
import os
import allel
import numpy as np

vcf_file = sys.argv[1]
sample_population_file = sys.argv[2]
window_size =  sys.argv[3]
window_slide = sys.argv[4]
output_file = sys.argv[5]

sample_population_map = {}
with open(sample_population_file, 'r') as pop_file:
    for line in pop_file:
        line = line.strip()
        values = line.split('\t')
        if len(values) == 2:
            sample_name, population = values
            #sample_name, population = line.strip().split('\t')
            sample_population_map[sample_name] = int(population)
        else:
            print(f"Ignoring line: {line}")

callset = allel.read_vcf(vcf_file, fields=['samples', 'variants/POS', 'calldata/GT'])

with open(output_file, 'w') as outfile:
    outfile.write("Start_Pos\tEnd_Pos\tAverage_FST\tWeighted_FST\n")
    for start_pos in range(1, callset['variants/POS'][-1], window_slide):
        end_pos = start_pos + window_size
        loc_region = callset['variants/POS'][(callset['variants/POS'] >= start_pos) & (callset['variants/POS'] < end_pos)]
        loc_index = allel.SortedIndex(loc_region)
        sub_vcf = allel.subset_vcf(callset, loc_index, fields=['calldata/GT'])
        population_gts = {}
        for pop_num in range(1, 101):  # Assuming populations are numbered from 1 to 100
            population_gts[pop_num] = []
        for sample_index, sample_name in enumerate(sub_vcf['samples']):
            population_num = sample_population_map[sample_name]
            population_gts[population_num].append(sub_vcf['calldata/GT'][:, sample_index])
        for pop_num in range(1, 101):
            population_gts[pop_num] = np.array(population_gts[pop_num])
        population_af = {}
        for pop_num in range(1, 101):
            population_af[pop_num] = allel.mean_allele_frequency(population_gts[pop_num])
        avg_fsts = []
        for i in range(1, 101):
            for j in range(i + 1, 101):
                fst_values = allel.weir_cockerham_fst(population_gts[i], population_gts[j])
                avg_fst = np.nanmean(fst_values)
                avg_fsts.append(avg_fst)
        weighted_fsts = []
        for i in range(1, 101):
            for j in range(i + 1, 101):
                n_samples_pop_i = population_gts[i].shape[1]
                n_samples_pop_j = population_gts[j].shape[1]
                total_samples = n_samples_pop_i + n_samples_pop_j
                weighted_fst = (n_samples_pop_i * n_samples_pop_j) / (total_samples - 1) * (1 - np.sum(population_af[i] * population_af[j]))
                weighted_fsts.append(weighted_fst)
        outfile.write(f"{start_pos}\t{end_pos}\t{np.nanmean(avg_fsts)}\t{','.join(map(str, weighted_fsts))}\n")

