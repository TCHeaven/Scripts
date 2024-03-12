import sys

fasta_file_path = sys.argv[1]
output_path = sys.argv[2]

def process_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    max_length = {}
    max_sequence = {}
    current_identifier = ""
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            current_identifier = line.split('|')[0][1:]
        else:
            if current_identifier not in max_length:
                max_length[current_identifier] = 0
            sequence_length = len(line)
            if sequence_length > max_length[current_identifier]:
                max_length[current_identifier] = sequence_length
                max_sequence[current_identifier] = line
    with open(output_path, 'w') as output_file:
        for identifier in max_length:
            output_file.write(f">{identifier}\n")
            output_file.write(f"{max_sequence[identifier]}\n")

process_fasta(fasta_file_path)

