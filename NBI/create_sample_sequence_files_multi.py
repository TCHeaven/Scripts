import sys
import os
import re
from multiprocessing import Pool

gene_info_file = sys.argv[1]
vcf_file = sys.argv[2]
reference_fasta = sys.argv[3]
output_directory = sys.argv[4]

print(gene_info_file)
print(vcf_file)
print(reference_fasta)
print(output_directory)

def parse_gene_info(gene_info_file):
    # This function parses the gene information file and returns a dictionary with gene information.
    gene_info = {}
    with open(gene_info_file, 'r') as file:
        for line in file:
            # Regular expression pattern to extract relevant information from each line
            match = re.match(r'Chromosome: (\S+)\s+Gene Name: (\S+)\s+Start Position: (\d+)\s+Stop Position: (\d+)', line)
            if match:
                chromosome = match.group(1)
                gene_name = match.group(2)
                start_position = int(match.group(3))
                stop_position = int(match.group(4))
                gene_info[gene_name] = {'chromosome': chromosome, 'start': start_position, 'stop': stop_position}
    return gene_info

def load_gene_sequence(reference_fasta, gene_name):
    # This function loads the gene sequence from the reference fasta file based on the gene name.
    fasta_file = f'{reference_fasta}'
    gene_sequence = ""
    with open(fasta_file, 'r') as file:
        lines = file.readlines()
        for i in range(len(lines)):
            line = lines[i].strip()
            if line.startswith('>') and gene_name in line:
                # Found the header line for the gene
                sequence_lines = []
                for j in range(i + 1, len(lines)):
                    next_line = lines[j].strip()
                    if next_line.startswith('>'):
                        # Reached the next gene, break the loop
                        break
                    sequence_lines.append(next_line)
                gene_sequence = "".join(sequence_lines)
                break
    return gene_sequence

def process_vcf_chunk(samples_chunk, gene_info, reference_fasta, output_directory):
    for gene_name, info in gene_info.items():
        gene_fasta = load_gene_sequence(reference_fasta, gene_name)
        # Create reference gene FASTA file
        reference_gene_file = f'{output_directory}/hom_{gene_name}_REFERENCE.fasta'
        with open(reference_gene_file, 'w') as reference_file:
            reference_file.write(f'>{gene_name}\n{gene_fasta}\n')
        # Process each sample in the chunk
        for sample_data in samples_chunk:
            sample_name, sample_lines = sample_data
            for line in sample_lines:
                if not line.startswith('#'):
                    fields = line.split('\t')
                    chromosome = fields[0]
                    position = int(fields[1])
                    ref_base = fields[3]
                    alt_bases = fields[4].split(',')
                    sample_values = fields[9:]
                    if info['chromosome'] == chromosome and info['start'] <= position <= info['stop']:
                        gene_start_position = info['start']
                        gene_relative_position = position - gene_start_position
                        print("gene_relative_position:", gene_relative_position)
                        print("length of gene_fasta:", len(gene_fasta))
                        gene_ref_base = gene_fasta[gene_relative_position]
                        if gene_ref_base != ref_base:
                            print(f"Mismatch in gene {gene_name} at position {position}, start {gene_start_position}, relative {gene_relative_position}: Expected {ref_base}, Found {gene_ref_base}")
                        # Create mutant gene FASTA files for each sample
                    for sample_name, sample_value in zip(sample_names, sample_values):
                        print(f"Sample: {sample_name}")
                        if sample_value.split(':')[0].split('/')[0] == '.':
                            print(f"SNP at {position} is missing in this sample")
                            alt_base = "N"
                        else:
                            allele_index = int(sample_value.split(':')[0].split('/')[0])
                            print(f"Allele index at {position}: {allele_index}")
                            if allele_index == 0:
                                    alt_base = ref_base
                            else:
                                alt_base = alt_bases[allele_index - 1] #always 0?
                        mutant_gene_fasta = gene_fasta[:gene_relative_position] + alt_base + gene_fasta[gene_relative_position + 1:]
                        mutant_gene_file = f'{output_directory}/hom_{gene_name}_{sample_name}_SAMPLE.fasta'
                        # Check if the mutant gene FASTA file already exists
                        if os.path.exists(mutant_gene_file):
                            # Read the existing file and update the SNP position with the alternate allele
                            with open(mutant_gene_file, 'r') as existing_file:
                                existing_header = existing_file.readline().strip()
                                existing_sequence = existing_file.readline().strip()
                            existing_sequence_list = list(existing_sequence)
                            existing_sequence_list[gene_relative_position] = alt_base
                            updated_sequence = "".join(existing_sequence_list)
                            # Overwrite the existing file with the updated sequence
                            with open(mutant_gene_file, 'w') as updated_file:
                                updated_file.write(f'{existing_header}\n{updated_sequence}\n')
                        else:
                            # Create a new mutant gene FASTA file
                            with open(mutant_gene_file, 'w') as mutant_file:
                                mutant_file.write(f'>{gene_name}_{sample_name}\n{mutant_gene_fasta}\n')

if __name__ == "__main__":
    # Parse gene information from the file
    gene_info = parse_gene_info(gene_info_file)
    # Split the VCF file into chunks based on samples to process concurrently
    num_processes = 4  # Number of available CPU cores
    samples_per_chunk = 10  # Adjust this value based on the number of samples you want to process per chunk
    chunks = []
    with open(vcf_file, 'r') as vcf:
        header_line = vcf.readline()
        sample_names = header_line.strip().split('\t')[9:]
        while True:
            try:
                line = next(vcf)
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    sample_info = fields[0:9]
                    sample_lines = [line]
                    for _ in range(samples_per_chunk - 1):
                        next_line = next(vcf)
                        if not next_line.startswith('#'):
                            sample_lines.append(next_line)
                    chunks.append([(sample_name, sample_info, sample_lines) for sample_name, sample_info in zip(sample_names, zip(*[sample_info]))])
            except StopIteration:
                break    # Create a pool of worker processes
    pool = Pool(processes=num_processes)
    # Process the VCF chunks using the pool of worker processes
    pool.starmap(process_vcf_chunk, [(chunk, gene_info, reference_fasta, output_directory) for chunk in chunks])
    # Close the pool of worker processes
    pool.close()
    pool.join()

