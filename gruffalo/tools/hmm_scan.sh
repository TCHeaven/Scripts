#!/bin/bash
#SBATCH --partition=himem
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH --time=240:00:00


CurPath=$PWD
HmmFile=$1
ProtFile=$2
Prefix=$3
OutDir=$4


WorkDir=$TMPDIR/hmmscan

mkdir -p $WorkDir
cd $WorkDir
cp $CurPath/$ProtFile proteins.fa
# cp $CurPath/$HmmFile hmm.txt

hmmscan --cpu 8 --domtblout $Prefix.out.dm $CurPath/$HmmFile proteins.fa > $Prefix.out
# hmmscan --cpu 16 --domtblout $Prefix.out.dm hmm.txt proteins.fa > $Prefix.out

mkdir -p $CurPath/$OutDir
cp $Prefix.out.dm $CurPath/$OutDir/.
cp $Prefix.out $CurPath/$OutDir/.


