import sys

input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

try:
    with open(input_file_path, 'r') as file:
        lines = file.readlines()
except FileNotFoundError:
    print(f"Error: File '{input_file_path}' not found.")
    sys.exit(1)

#Remove everything after ":" in columns 1 and 4
processed_lines = []
for line in lines:
    columns = line.strip().split('\t')
    columns[0] = columns[0].split(':')[0]  # Remove everything after ":" in column 1
    columns[3] = columns[3].split(':')[0]  # Remove everything after ":" in column 4
    processed_lines.append('\t'.join(columns))

with open(output_file_path, 'w') as file:
    file.write('\n'.join(processed_lines))

print(f"Processed data has been saved to '{output_file_path}'.")

