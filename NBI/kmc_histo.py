import sys
import struct

kmc_file = sys.argv[1]
output_file = sys.argv[2]

occurrences = [0] * 100001

with open(kmc_file, 'r') as file:
    for line in file:
        _, number = line.strip().split()
        number = int(number)
        if 1 <= number <= 100000:
            occurrences[number] += 1

with open(output_file, 'w') as file:
    for i, count in enumerate(occurrences[1:], start=1):
        file.write(f'{i}\t{count}\n')
