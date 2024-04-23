import sys

def convert_to_single_line(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # If the line is a header, write it directly to the output file
                outfile.write(line)
            else:
                # If the line is sequence data, strip newline characters and write to the output file
                outfile.write(line.strip())


if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python script.py input.fasta output.fasta")
        sys.exit(1)

    # Extract input and output file paths from command-line arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Call the conversion function with the provided file paths
    convert_to_single_line(input_file, output_file)
