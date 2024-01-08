import sys
import os
import re

input_file_list = sys.argv[1]

def filter_common_lines(input_file_list, output_prefix='non-zero_'):
    def create_id_set_and_filter(lines):
        id_set = set()
        filtered_lines = []
        for line in lines:
            if line.startswith("#") or line.startswith("chr"):
                continue            
            columns = line.strip().split('\t')
            line_id = (columns[0], int(columns[1]), columns[2])
            try:
                ratio = float(columns[4])
            except ValueError:
                continue  # Skip lines with non-numeric or invalid values
            if ratio != 0.000:
                id_set.add(line_id)
                filtered_lines.append(line)
        return id_set, filtered_lines
    with open(input_file_list, 'r') as infile:
        input_files = [line.strip() for line in infile.readlines()]
    id_sets, filtered_lines_list = zip(*[create_id_set_and_filter(open(input_file, 'r').readlines()) for input_file in input_files])
    common_ids = set.intersection(*id_sets)
    for input_file, filtered_lines in zip(input_files, filtered_lines_list):
        output_file = '/'.join([output_prefix + input_file.split('/')[-1]])
        with open(output_file, 'w') as out:
            for line in filtered_lines:
                out.write(line)



filter_common_lines(input_file_list)
