from Bio import SeqIO
import sys
import os

signalp_file = sys.argv[1]
input_fasta = sys.argv[2]
output = os.path.splitext(input_fasta)[0]

output_fasta = f"{output}-cleaved.faa"
output_directory = os.path.splitext(output_fasta)[0]

def load_signal_peptide_info(signalp_file):
    cleavage_sites = {}
    with open(signalp_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            seq_id = parts[0]
            cleavage_site = int(parts[1])
            cleavage_sites[seq_id] = cleavage_site
    return cleavage_sites

def remove_signal_peptide_fasta(input_fasta, output_fasta, cleavage_sites):
    with open(output_fasta, 'w') as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq_id = record.id.split('.')[0]  # To match with IDs in cleavage_sites
            if seq_id in cleavage_sites:
                cleavage_site = cleavage_sites[seq_id]
                mature_sequence = record.seq[cleavage_site:]
                record.seq = mature_sequence
            SeqIO.write(record, output_handle, "fasta")

cleavage_sites = load_signal_peptide_info(signalp_file)
remove_signal_peptide_fasta(input_fasta, output_fasta, cleavage_sites)

def split_fasta(output_fasta, output_directory):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    for record in SeqIO.parse(output_fasta, "fasta"):
        output_file = os.path.join(output_directory, f"{record.id}.fasta")
        with open(output_file, 'w') as output_handle:
            SeqIO.write(record, output_handle, "fasta")
        print(f"Written {record.id} to {output_file}")

split_fasta(output_fasta, output_directory)
