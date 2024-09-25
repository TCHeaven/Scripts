#!/usr/bin/python

import sys
filename = sys.argv[1]
with open(filename, 'r') as file:
    lines = file.readlines()

five_prime_values = []

for line in lines[1:]:
    columns = line.strip().split()
    name = columns[0]
    if columns[1] == 'NA':
        continue
    five_prime = int(columns[1])
    strand = columns[3]
    five_prime_neighbour_strand = columns[4]
    if five_prime_neighbour_strand == 'NA':
        continue
    if five_prime != 99999 and strand != five_prime_neighbour_strand:
        five_prime_values.append(five_prime)

if five_prime_values:
    average_five_prime = sum(five_prime_values) / float(len(five_prime_values))
    print("Average different five_prime value:", average_five_prime)
else:
    print("No valid rows found for averaging.")

with open(filename, 'r') as file:
    lines = file.readlines()

five_prime_values = []

for line in lines[1:]:
    columns = line.strip().split()
    name = columns[0]
    if columns[1] == 'NA':
        continue
    five_prime = int(columns[1])
    strand = columns[3]
    five_prime_neighbour_strand = columns[4]
    if five_prime_neighbour_strand == 'NA':
        continue
    if five_prime != 99999 and strand == five_prime_neighbour_strand:
        five_prime_values.append(five_prime)

# Calculate the average of the filtered values
if five_prime_values:
    average_five_prime = sum(five_prime_values) / float(len(five_prime_values))
    print("Average same five_prime value:", average_five_prime)
else:
    print("No valid rows found for averaging.")

with open(filename, 'r') as file:
    lines = file.readlines()

three_prime_values = []

for line in lines[1:]:
    columns = line.strip().split()
    name = columns[0]
    if columns[2] == 'NA':
        continue
    three_prime = int(columns[2])
    strand = columns[3]
    three_prime_neighbour_strand = columns[5]
    if three_prime_neighbour_strand == 'NA':
        continue
    if three_prime != 99999 and strand != three_prime_neighbour_strand:
        three_prime_values.append(three_prime)

if three_prime_values:
    average_three_prime = sum(three_prime_values) / float(len(three_prime_values))
    print("Average different three_prime value:", average_three_prime)
else:
    print("No valid rows found for averaging.")

with open(filename, 'r') as file:
    lines = file.readlines()

three_prime_values = []

for line in lines[1:]:
    columns = line.strip().split()
    name = columns[0]
    if columns[2] == 'NA':
        continue
    three_prime = int(columns[2])
    strand = columns[3]
    three_prime_neighbour_strand = columns[5]
    if three_prime_neighbour_strand == 'NA':
        continue
    if three_prime != 99999 and strand == three_prime_neighbour_strand:
        three_prime_values.append(three_prime)

if three_prime_values:
    average_three_prime = sum(three_prime_values) / float(len(three_prime_values))
    print("Average same three_prime value:", average_three_prime)
else:
    print("No valid rows found for averaging.")
