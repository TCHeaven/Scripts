#!/usr/bin/env bash

# This script combines minimap2, miniasm and minipolish into a single command
# for easier execution.

# It takes two positional arguments:
#  1) a long read file
#  2) the number of threads to use


# Create temporary intermediate files.
overlaps=$(mktemp)".paf"
unpolished_assembly=$(mktemp)".gfa"

# Find read overlaps with minimap2.
minimap2 -x ava-ont -t "$2" "$1" "$1" > "$overlaps"

# Run miniasm to make an unpolished assembly.
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven \
/data/users/theaven/miniasm_0.3--h577a1d6_5 \
miniasm -f "$1" "$overlaps" > "$unpolished_assembly"

# Polish the assembly with minipolish, outputting the result to stdout.
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven \
/data/users/theaven/minipolish_0.2.1--pyhdfd78af_0 \
minipolish --threads "$2" "$1" "$unpolished_assembly"

# Clean up.
rm "$overlaps" "$unpolished_assembly"
