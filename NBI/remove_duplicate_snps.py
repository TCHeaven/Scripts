import sys
import os

vcf_file = sys.argv[1]

# Read the VCF file
with open(vcf_file, 'r') as input_file:
    lines = input_file.readlines()

# Create a list to store unique SNP lines
unique_snps = []

# Create a set to store unique SNP identifiers
unique_ids = set()

# Iterate over each line in the VCF file
for line in lines:
    if line.startswith('#'):  # Skip header lines
        unique_snps.append(line)
    else:
        # Extract chromosome, position, reference allele, and alternate allele
        fields = line.split('\t')
        chromosome = fields[0]
        position = fields[1]
        ref_allele = fields[3]
        alt_allele = fields[4]
        snp_id = f"{chromosome}_{position}_{ref_allele}_{alt_allele}"

        if snp_id not in unique_ids:
            unique_snps.append(line)
            unique_ids.add(snp_id)

# Write the unique SNP lines to a new VCF file
output_filename = f'dedup_{vcf_file}'
with open(output_filename, 'w') as output_file:
    output_file.writelines(unique_snps)

print(f"Unique SNPs written to {output_filename}")

