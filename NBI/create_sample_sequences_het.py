import sys
import os
import gzip
from collections import defaultdict

vcf_file = sys.argv[1]
reference_fasta = sys.argv[2]
output_dir = sys.argv[3]

header_line = ''
sample_snps = defaultdict(list)

with gzip.open(vcf_file, 'rt') as file:
    for line in file:
        if line.startswith('#CHROM'):
            header_line = line.strip()
            sample_names = header_line.split('\t')[9:]
            break
    for line in file:
        if not line.startswith('#'):
            fields = line.strip().split('\t')
            chromosome = fields[0]
            position = int(fields[1])
            ref_base = fields[3]
            alt_bases = fields[4]
            sample_values = fields[9:]
            for i, sample_value in enumerate(sample_values):
                sample_name = sample_names[i]
                genotype = sample_value.split(':')[0]
                genotype_lookup = {
                    '0/0': ref_base,
                    './.': ref_base,
                    '0/1': alt_bases[0],
                    '1/0': alt_bases[0],
                    '1/1': alt_bases[0]
                }
                allele = genotype_lookup.get(genotype)
                if allele is None:
                    print(f'ERROR with VCF format: {sample_value}')
                else:
                    sample_snps[sample_name].append((chromosome, position, allele, ref_base))

chromosome_sequences = {}
current_chromosome = ''

with open(reference_fasta, 'r') as fasta_file:
    sequence_lines = []
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            if current_chromosome:
                chromosome_sequences[current_chromosome] = ''.join(sequence_lines)
                sequence_lines = []
            current_chromosome = line[1:]
        else:
            sequence_lines.append(line)
    chromosome_sequences[current_chromosome] = ''.join(sequence_lines)

for sample_name in sample_snps:
    print(f"{sample_name}")
    fasta_filename = f'{output_dir}/{sample_name}_het.fasta'
    with open(fasta_filename, 'w') as fasta_file:
        sequence = chromosome_sequences[sample_snps[sample_name][0][0]]
        current_chromosome = None
        for chromosome, position, allele, ref_base in sample_snps[sample_name]:
            if chromosome != current_chromosome:
                if current_chromosome is not None:
                    fasta_file.write(f'>{current_chromosome}\n{sequence}\n')
                current_chromosome = chromosome
                sequence = chromosome_sequences[current_chromosome]
            if ref_base == sequence[position - 1]:
                #print(sequence[position - 3])
                #print(sequence[position - 2])
                #print(sequence[position - 1])
                #print(sequence[position])
                #print(sequence[position + 1])
                #print(sequence[position + 2])
                updated_sequence = sequence[:position - 1] + allele + sequence[position:]
                sequence = updated_sequence
                if allele == sequence[position - 1]:
                    print(f"Replaced")
                    print(f"{allele}")
                    #print(sequence[position - 3])
                    #print(sequence[position - 2])
                    #print(sequence[position - 1])
                    #print(sequence[position])
                    #print(sequence[position + 1])
                    #print(sequence[position + 2])
                else:
                    print(f"Not Replaced")
            else:
                actual = sequence[position]
                print(f'Warning: Reference base {ref_base} in VCF does not match reference fasta {actual} at {chromosome}:{position}')
        fasta_file.write(f'>{current_chromosome}\n{sequence}\n')

print("Fasta files generated successfully.")
#import sys
#import os
#import gzip

#vcf_file = sys.argv[1]
#reference_fasta = sys.argv[2]
#output_dir = sys.argv[3]

#header_line = ''
#sample_names = []
#sample_snps = {}

#with gzip.open(vcf_file, 'rt') as file:
#    for line in file:
#        if line.startswith('#CHROM'):
#            header_line = line.strip()
#            sample_names = header_line.split('\t')[9:]
#            for sample_name in sample_names:
#                sample_snps[sample_name] = []
#            break
#    for line in file:
#        if not line.startswith('#'):
#            fields = line.strip().split('\t')
#            chromosome = fields[0]
#            position = int(fields[1])
#            ref_base = fields[3]
#            alt_bases = fields[4]
#            sample_values = fields[9:]    
#            for i, sample_value in enumerate(sample_values):
#                sample_name = sample_names[i]
#                genotype = sample_value.split(':')[0]
#                if genotype == '0/0':
#                    allele = ref_base
#                elif genotype == './.':
#                    allele = ref_base
#                elif genotype == '0/1' or genotype == '1/0':
#                    allele = ref_base
#                elif genotype == '1/1':
#                    allele = alt_bases[0] 
#                else:
#                    print(f'ERROR with vcf format')
#                    print(sample_value)
#                sample_snps[sample_name].append((chromosome, position, allele, ref_base))
#
#chromosome_sequences = {}
#current_chromosome = ''
#with open(reference_fasta, 'r') as fasta_file:
#    sequence_lines = []
#    for line in fasta_file:
#        line = line.strip()
#        if line.startswith('>'):
#            if current_chromosome:
#                chromosome_sequences[current_chromosome] = ''.join(sequence_lines)
#                sequence_lines = []
#            current_chromosome = line[1:]
#        else:
#            sequence_lines.append(line)
#    chromosome_sequences[current_chromosome] = ''.join(sequence_lines)
#
#for sample_name in sample_snps:
#    print(f"{sample_name}")
#    sorted_snps = sorted(sample_snps[sample_name], key=lambda x: x[0])
#    fasta_filename = f'{output_dir}{sample_name}_hom.fasta'
#    with open(fasta_filename, 'w') as fasta_file:
#        sequence = chromosome_sequences[sample_snps[sample_name][0][0]]  
#        current_chromosome = None  
#        for chromosome, position, allele, ref_base in sample_snps[sample_name]:
#            if chromosome != current_chromosome:
#                if current_chromosome == None:
#                    print(f"Starting {sample_name}")
#                    current_chromosome = chromosome  # Update the current chromosome
#                    sequence = chromosome_sequences[current_chromosome]
#                    #print(f"Writing chromosome: {current_chromosome}")
#                    if ref_base == sequence[position - 1]:
#                        #print(f"1st")
#                        updated_sequence = sequence[:position - 2] + allele + sequence[position - 1:]
#                        sequence = updated_sequence
#                        #if allele == sequence[position - 1]:
#                            #print(f"Replaced")
#                        #else:
#                            #print(f"Not Replaced")
#                    else:
#                        actual = sequence[position]
#                        print(f'Warning: Reference base {ref_base} in VCF does not match reference fasta {actual} at {chromosome}:{position}')
#                else:
#                    fasta_file.write(f'>{current_chromosome}\n{sequence}\n')
#                    #print(f"Report")
#                    current_chromosome = chromosome  # Update the current chromosome
#                    sequence = chromosome_sequences[current_chromosome]
#                    #print(f"Writing chromosome: {current_chromosome}")
#                    if ref_base == sequence[position - 1]:
#                        #print(f"1st")
#                        updated_sequence = sequence[:position - 2] + allele + sequence[position - 1:]
#                        sequence = updated_sequence
#                        #if allele == sequence[position - 1]:
#                            #print(f"Replaced")
#                        #else:
#                            #print(f"Not Replaced")
#                    else:
#                        actual = sequence[position]
#                        print(f'Warning: Reference base {ref_base} in VCF does not match reference fasta {actual} at {chromosome}:{position}')
#            else:
#                if ref_base == sequence[position - 1]:
#                    #print(f"2nd")
#                    updated_sequence = sequence[:position - 2] + allele + sequence[position - 1:]
#                    sequence = updated_sequence
#                    #if allele == sequence[position - 1]:
#                        #print(f"Replaced")
#                    #else:
#                        #print(f"Not Replaced")
#                else:
#                    actual = sequence[position]
#                    print(f'Warning: Reference base {ref_base} in VCF does not match reference fasta {actual} at {chromosome}:{position}')
#        fasta_file.write(f'>{current_chromosome}\n{sequence}\n')
#        #print(f"Report")
#
#print("Fasta files generated successfully.")
