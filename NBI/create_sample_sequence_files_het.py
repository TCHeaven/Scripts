import sys
import os
import re

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

def process_vcf(vcf_file, gene_info, reference_fasta, output_directory):
    # This function processes the VCF file and creates sample gene FASTA files.
    with open(vcf_file, 'r') as file:
        header_line = ''
        sample_names = []
        for line in file:
            if line.startswith('#CHROM'):
                header_line = line.strip()
                sample_names = header_line.split('\t')[9:]
                break
        for line in file:
            if not line.startswith('#'):
                fields = line.split('\t')
                chromosome = fields[0]
                position = int(fields[1])
                ref_base = fields[3]
                alt_bases = fields[4].split(',')
                sample_values = fields[9:]
                for gene_name, info in gene_info.items():
                    if info['chromosome'] == chromosome and info['start'] <= position <= info['stop']:
                        gene_start_position = info['start']
                        gene_relative_position = position - gene_start_position
                        gene_fasta = load_gene_sequence(reference_fasta, gene_name)
                        gene_ref_base = gene_fasta[gene_relative_position]
                        if gene_ref_base != ref_base:
                            print(f"Mismatch in gene {gene_name} at position {position}: Expected {ref_base}, Found {gene_ref_base}")
                        # Create reference gene FASTA file
                        reference_gene_file = f'{output_directory}/{gene_name}_REFERENCE.fasta'
                        with open(reference_gene_file, 'w') as reference_file:
                            reference_file.write(f'>{gene_name}\n{gene_fasta}\n')
                        # Create mutant gene FASTA files for each sample
                        for sample_name, sample_value in zip(sample_names, sample_values):
                            allele_index = int(sample_value.split(':')[0].split('/')[1])
                            print(sample_name)
                            print(allele_index)
                            if allele_index == 0:
                                    alt_base = ref_base
                            else:
                                    alt_base = alt_bases[allele_index - 1]
                            mutant_gene_fasta = gene_fasta[:gene_relative_position] + alt_base + gene_fasta[gene_relative_position + 1:]
                            mutant_gene_file = f'{output_directory}/{gene_name}_{sample_name}_SAMPLE.fasta'
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


# Parse gene information from the file
gene_info = parse_gene_info(gene_info_file)

# Process the VCF file and create mutant gene fasta files for each sample
process_vcf(vcf_file, gene_info, reference_fasta, output_directory)

