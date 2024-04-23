from collections import defaultdict
import sys
import os

# Check if the correct number of arguments is provided
if len(sys.argv) != 2:
    print("Usage: python script.py input_file")
    sys.exit(1)

input_file = sys.argv[1]
output_file = os.path.splitext(input_file)[0] + "_mperc.csv"

# Load distance matrix
acc = defaultdict(list)
acc_order = []
with open(input_file) as inp:
    next(inp)
    for line in inp:
        A = line.strip().split()
        acc[A[0]] = A[1:]
        acc_order.append(A[0])

# Write the distance matrix in a desired format
with open(output_file, "w") as outf:
    outf.write(",".join(["ID"] + acc_order) + "\n")
    for a in acc_order:
        outf.write(",".join([a] + acc[a]) + "\n")

print("Done")
