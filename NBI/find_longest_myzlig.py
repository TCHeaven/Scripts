import sys
import os

fasta_file_path = sys.argv[1]
output_path = sys.argv[2]

def process_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    max_length_myz = 0
    max_sequence_myz = ""
    max_length_lig = 0
    max_sequence_lig = ""
    current_header = ""
    for line in lines:
        line = line.strip()
        if line.startswith('>MYZ') or line.startswith('>LIG'):
            current_header = line
        else:
            sequence_length = len(line)
            if current_header.startswith('>MYZ') and sequence_length > max_length_myz:
                max_length_myz = sequence_length
                max_sequence_myz = line
                max_header_myz = current_header
            elif current_header.startswith('>LIG') and sequence_length > max_length_lig:
                max_length_lig = sequence_length
                max_sequence_lig = line
                max_header_lig = current_header
    with open(output_path, 'w') as output_file:
        output_file.write(f"{max_header_myz}\n")
        output_file.write(f"{max_sequence_myz}\n")
        output_file.write(f"{max_header_lig}\n")
        output_file.write(f"{max_sequence_lig}\n")

process_fasta(fasta_file_path)
