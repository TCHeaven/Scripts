import os
import sys
import pandas as pd

input_file = sys.argv[1]
window_size = sys.argv[2]
output_file = sys.argv[3]

# Read the input data
data = pd.read_csv(input_file, sep='\t', header=None, names=["chromosome", "start", "end", "class", "score", "strand"])

# Initialize a dictionary to store TE densities for each chromosome
te_density = defaultdict(lambda: defaultdict(int))

# Process each TE annotation and add its length to the corresponding windows
for index, row in data.iterrows():
    chromosome = row['chromosome']
    start_window = row['start'] // window_size
    end_window = row['end'] // window_size
    for window in range(start_window, end_window + 1):
        window_start = window * window_size
        window_end = window_start + window_size
        # Calculate the overlap between the TE and the current window
        overlap_start = max(row['start'], window_start)
        overlap_end = min(row['end'], window_end)
        overlap_length = overlap_end - overlap_start
        if overlap_length > 0:
            te_density[chromosome][window] += overlap_length

# Write the results to the output file
with open(output_file, 'w') as out:
    for chromosome, windows in te_density.items():
        for window, total_te_length in sorted(windows.items()):
            window_start = window * window_size
            window_end = window_start + window_size
            density = total_te_length / window_size
            out.write(f"{chromosome}\t{window_start}\t{window_end}\t{density}\n")

print(f"Overall TE density data has been written to {output_file}")
