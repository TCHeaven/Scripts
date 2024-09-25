from Bio import Phylo
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

def collapse_single_leaf_nodes(tree):
    for clade in tree.find_clades(order='postorder'):
        if len(clade.clades) == 1 and clade.clades[0].is_terminal():
            leaf = clade.clades[0]
            clade.branch_length = (clade.branch_length or 0) + (leaf.branch_length or 0)
            clade.name = leaf.name
            clade.clades = []

def main():
    tree = Phylo.read(input_file, "newick")
    collapse_single_leaf_nodes(tree)
    Phylo.write(tree, output_file, "newick")

if __name__ == "__main__":
    main()
