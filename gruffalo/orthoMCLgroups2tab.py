#!/usr/bin/python3

import sys

class Cluster:
    def __init__(self, name):
        self.name = name
        self.members = {}

    def add(self, strain, member):
        self.members[strain] = member

all_proteins = {}
clusters = []
strains = {}

# Open the first file and process proteins
with open(sys.argv[1]) as file:
    for ln in file:
        if ln.startswith('>'):
            protein = ln[1:].rstrip()
            all_proteins[protein] = False

# Open the second file and process clusters
with open(sys.argv[2]) as file:
    for ln in file:
        prefix, proteins = ln.rstrip().split(':')
        proteins = proteins.split()

        c = Cluster(prefix)

        for p in proteins:
            strain, id = p.split("|")
            c.add(strain, id)
            strains[strain] = strain

            all_proteins[p] = True

        clusters.append(c)

# Add the singletons
i = 1
for p, already_counted in all_proteins.items():
    if already_counted == False:
        c = Cluster(f"single{i}")
        strain, id = p.split("|")
        c.add(strain, id)
        i += 1
        clusters.append(c)

# Output the names of clusters
first_cluster = clusters[0]
print(" ".join([f'"{c.name}"' for c in clusters]))

# Output the membership matrix
for b in sorted(strains):
    print(f'"{b}"', end=" ")

    for c in clusters:
        if b in c.members:
            print("1", end=" ")
        else:
            print("0", end=" ")
    print()
