import sys

input_fasta_file = sys.argv[1]

current_chr = ""
current_length = 0

with open(input_fasta_file, "r") as fasta_file:
    for line in fasta_file:
        if line.startswith(">"):
            # If it's a header line, output the previous chromosome info
            if current_chr:
                print(f"{current_chr}\t{current_length}")
            # Extract chromosome name from the header
            current_chr = line.strip().lstrip(">")
            current_length = 0
        else:
            # Count the number of base pairs in the sequence
            current_length += len(line.strip())

# Output the last chromosome info
if current_chr:
    print(f"{current_chr}\t{current_length}")

