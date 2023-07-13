import sys
import os
import gzip
from collections import defaultdict
from multiprocessing import Pool

# Read command-line arguments
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
                    '0/1': ref_base,
                    '1/0': ref_base,
                    '1/1': alt_bases[0]
                }
                allele = genotype_lookup.get(genotype)
                if allele is None:
                    print(f'ERROR with VCF format: {sample_value}')
                else:
                    sample_snps[sample_name].append((chromosome, position, allele, ref_base))

print("VCF loaded")

chromosome_sequences = {}
current_chromosome = ''

# Load chromosome sequences into memory only once
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

print("Fasta loaded")

def process_sample_snps(sample_name):
    print(f"Processing sample: {sample_name}")
    fasta_filename = f'{output_dir}/{sample_name}_hom2.fasta'
    updated_sequences = {}
    sequence = chromosome_sequences[sample_snps[sample_name][0][0]]
    current_chromosome = None
    for chromosome, position, allele, ref_base in sample_snps[sample_name]:
        if chromosome != current_chromosome:
            if current_chromosome is not None:
                updated_sequences[current_chromosome] = sequence
            current_chromosome = chromosome
            sequence = chromosome_sequences[current_chromosome]
        if ref_base == sequence[position - 1]:
            updated_sequence = sequence[:position - 1] + allele + sequence[position:]
            sequence = updated_sequence
    updated_sequences[current_chromosome] = sequence
    with open(fasta_filename, 'w') as fasta_file:
        for chromosome, sequence in updated_sequences.items():
            fasta_file.write(f'>{chromosome}\n{sequence}\n')
    print(f"Finished processing sample: {sample_name}")

print(f"CPUs: {os.cpu_count()}")

if __name__ == '__main__':
    # Use multiprocessing to process samples concurrently
    with Pool(processes=8) as pool:
        pool.map(process_sample_snps, sample_snps.keys())

    print("Fasta files generated successfully.")
