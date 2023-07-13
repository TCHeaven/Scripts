import sys
import os
import pysam
from pyfaidx import Fasta

# Provide the file paths for the VCF file, gene prediction file, and output directory
vcf_file = sys.argv[1]
gff3_file = sys.argv[2]
output_directory = sys.argv[3]


def extract_gene_snps(vcf_file, gff3_file, output_directory):
    # Parse the GFF3 file to extract gene coordinates
    gene_coordinates = {}
    with open(gff3_file, "r") as gff_file:
        for line in gff_file:
            if line.startswith("#"):
                continue
            gene_info = line.strip().split("\t")
            if gene_info[2] != "CDS":
                continue
            gene_id = gene_info[8].split(";")[0].replace("ID=", "").split(".")[0:2]
            gene_chrom = gene_info[0]  # Store the chromosome name from the GFF3 file
            gene_start = int(gene_info[3])
            gene_end = int(gene_info[4])
            gene_coordinates[gene_id] = (gene_chrom, gene_start, gene_end)
    # Parse the VCF file and write gene-specific SNP data to separate VCF files
    with open(vcf_file, "r") as vcf:
        gene_vcf_files = {}
        closed_files = []  # List to store closed gene-specific VCF files
        header_line = ""  # Variable to store the header line
        num_open_files = 0  # Counter for the number of open files
        for line in vcf:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    header_line = line.strip()  # Store the header line
                continue
            vcf_info = line.strip().split("\t")
            chrom = vcf_info[0]
            pos = int(vcf_info[1])
            # Check if SNP position falls within any gene
            for gene_id, (gene_chrom, gene_start, gene_end) in gene_coordinates.items():
                if chrom == gene_chrom and gene_start <= pos <= gene_end:
                    # Create a VCF file for the gene if it doesn't exist or has been closed
                    if gene_id not in gene_vcf_files or gene_id in closed_files:
                        if num_open_files >= 1000:
                            # Close the earliest opened file
                            oldest_file = next(iter(gene_vcf_files))
                            gene_vcf_files[oldest_file].close()
                            del gene_vcf_files[oldest_file]
                            closed_files.append(oldest_file)
                            num_open_files -= 1
                        if gene_id in closed_files:
                            # Reopen the closed gene-specific VCF file in append mode
                            gene_vcf_files[gene_id] = open(f"{output_directory}/{gene_id}_snps.vcf", "a")
                            closed_files.remove(gene_id)  # Remove the file from the closed files list
                        else:
                            gene_vcf_files[gene_id] = open(f"{output_directory}/{gene_id}_snps.vcf", "w")
                            gene_vcf_files[gene_id].write("##fileformat=VCFv4.3\n")
                            gene_vcf_files[gene_id].write(header_line + "\n")  # Write the header line
                            num_open_files += 1
                            print(f"Created file: {gene_id}_snps.vcf")
                     # Write the SNP information to the gene-specific VCF file
                    gene_vcf_files[gene_id].write("\t".join(vcf_info) + "\n")
                    print(f"Added SNP to {gene_id}_snps.vcf")                   
    # Close the remaining gene-specific VCF files
    for gene_id, vcf_file in gene_vcf_files.items():
        vcf_file.close()

# Call the function to extract gene-specific SNP data and create individual VCF files
extract_gene_snps(vcf_file, gff3_file, output_directory)
