from Bio import SeqIO
import sys

def get_longest_splice_variant(fasta_file, output_file):
    # Dictionary to store sequences grouped by gene
    gene_sequences = {}
    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Extract gene ID from the record ID
        gene_id = record.id.split(".t")[0]
        # Update gene_sequences dictionary
        if gene_id not in gene_sequences:
            gene_sequences[gene_id] = [record]
        else:
            gene_sequences[gene_id].append(record)
    # Dictionary to store the longest splice variant for each gene
    longest_splice_variants = {}
    # Iterate over gene_sequences
    for gene_id, sequences in gene_sequences.items():
        # Find the longest splice variant
        longest_sequence = max(sequences, key=len)
        longest_splice_variants[gene_id] = longest_sequence
    # Write longest splice variants to a new FASTA file
    with open(output_file, "w") as output_handle:
        SeqIO.write(longest_splice_variants.values(), output_handle, "fasta")

# Example usage:
fasta_file = sys.argv[1]
output_file = sys.argv[2]

get_longest_splice_variant(fasta_file, output_file)

