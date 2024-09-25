import os
import sys
import pandas as pd
from collections import defaultdict

if len(sys.argv) != 7:
  print("Usage: python script.py <input_file> <window_size> <feature> <karyotype> <output_dir> <output_prefix>")
  sys.exit(1)

input_file = sys.argv[1]
window_size = int(sys.argv[2]) 
feature = sys.argv[3]
karyotype = sys.argv[4]
output_dir = sys.argv[5]
output_prefix = sys.argv[6]

if not os.path.exists(output_dir):
  os.makedirs(output_dir)

def load_chromosome_lengths(karyotype):
  lengths = {}
  with open(karyotype, 'r') as f:
    for line in f:
      parts = line.strip().split()
      chromosome = parts[2] 
      chr_length = int(parts[5]) 
      lengths[chromosome] = chr_length
  return lengths

chromosome_lengths = load_chromosome_lengths(karyotype)

data = pd.read_csv(input_file, sep='\t', header=None, names=["chromosome", "tool", "feature", "start", "end", "score", "strand", "null", "description"])

density = defaultdict(lambda: defaultdict(int))

if feature == "all":
  for index, row in data.iterrows():
    chromosome = row['chromosome']
    start_window = row['start'] // window_size
    end_window = row['end'] // window_size
    for window in range(start_window, end_window + 1):
      window_start = window * window_size
      window_end = window_start + window_size
      overlap_start = max(row['start'], window_start)
      overlap_end = min(row['end'], window_end)
      overlap_length = overlap_end - overlap_start
      if overlap_length > 0:
        density[chromosome][window] += overlap_length
else:
  for index, row in data.iterrows():
    chromosome = row['chromosome']
    start_window = row['start'] // window_size
    end_window = row['end'] // window_size
    if row['feature'] == feature:
      for window in range(start_window, end_window + 1):
        window_start = window * window_size
        window_end = window_start + window_size
        overlap_start = max(row['start'], window_start)
        overlap_end = min(row['end'], window_end)
        overlap_length = overlap_end - overlap_start
        if overlap_length > 0:
          if feature in row['feature']:
            density[chromosome][window] += overlap_length

output_density = os.path.join(output_dir, f"{output_prefix}_density.tsv")
with open(output_density, 'w') as out:
  for chromosome in sorted(density.keys()):
    windows = density[chromosome]
    chromosome_length = chromosome_lengths.get(chromosome)
    for window in range((max(windows.keys()) + 1) if windows else 0):
      window_start = window * window_size
      window_end = window_start + window_size
      if window_end > chromosome_length:
        window_end = chromosome_length
      total_length = windows.get(window, 0)
      density_no = total_length / (window_end - window_start) if (window_end - window_start) > 0 else 0
      out.write(f"{chromosome}\t{window_start}\t{window_end}\t{density_no}\n")

print(f"Density data has been written to {output_density}")

count = defaultdict(lambda: defaultdict(int))

if feature == "all":
  for index, row in data.iterrows():
    chromosome = row['chromosome']
    start_window = row['start'] // window_size
    end_window = row['end'] // window_size
    for window in range(start_window, end_window + 1):
      window_start = window * window_size
      window_end = window_start + window_size
      overlap_start = max(row['start'], window_start)
      overlap_end = min(row['end'], window_end)
      overlap_length = overlap_end - overlap_start
      if overlap_length > 0:
        count[chromosome][window] += 1
else:
  for index, row in data.iterrows():
    chromosome = row['chromosome']
    start_window = row['start'] // window_size
    end_window = row['end'] // window_size
    if row['feature'] == feature:  
      for window in range(start_window, end_window + 1):
        window_start = window * window_size
        window_end = window_start + window_size
        overlap_start = max(row['start'], window_start)
        overlap_end = min(row['end'], window_end)
        overlap_length = overlap_end - overlap_start
        if overlap_length > 0:
          if feature in row['feature']:
            count[chromosome][window] += 1

output_count = os.path.join(output_dir, f"{output_prefix}_count.tsv")
with open(output_count, 'w') as out:
  for chromosome in sorted(count.keys()):
    windows = count[chromosome]
    chromosome_length = chromosome_lengths.get(chromosome)
    for window in range((max(windows.keys()) + 1) if windows else 0):
      window_start = window * window_size
      window_end = window_start + window_size
      if window_end > chromosome_length:
        window_end = chromosome_length
      total_count = windows.get(window, 0)
      out.write(f"{chromosome}\t{window_start}\t{window_end}\t{total_count}\n")

print(f"Number data has been written to {output_count}")
