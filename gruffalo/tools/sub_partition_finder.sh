#!/usr/bin/env bash
#SBATCH -J Partion_finder
#SBATCH --partition=long
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task 16


#### Submit a Partition Finder job in a folder specified 

input=$1

python=/mnt/shared/scratch/theaven/apps/conda/envs/predector2.7/bin/python2.7 
partition_finder=/home/theaven/scratch/apps/phylogeny/partitionfinder-2.1.1/PartitionFinderProtein.py
 
 $python $partition_finder $input
