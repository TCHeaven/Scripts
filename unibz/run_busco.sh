#!/bin/bash
#SBATCH -o /home/clusterusers/theaven/slurm_records/slurm.%j.out
#SBATCH -e /home/clusterusers/theaven/slurm_records/slurm.%j.out
#SBATCH --mem 32G
#SBATCH --nodes=1
#SBATCH --cpus-per-task 8
#SBATCH --account=shame
#SBATCH --partition=bioagri
#SBATCH --time=1-00:00:00

CurPath=$PWD
WorkDir="${TMPDIR:-/tmp}/${SLURM_JOB_ID}"
Genome="${1:?ERROR: Missing Genome}"
Database="${2:?ERROR: Missing Database}"
OutDir="${3:?ERROR: Missing Output Directory}"
OutFile="${4:?ERROR: Missing Output File}"
cpu="${SLURM_CPUS_PER_TASK:-1}"

echo CurPth:
echo $CurPath
echo WorkDir:
echo $WorkDir
echo OutDir:
echo $OutDir
echo OutFile:
echo $OutFile
echo Genome:
echo $Genome
echo Database:
echo $Database
echo _
echo _

mkdir -p $WorkDir
cp $Genome $WorkDir/genome.fa

cd $WorkDir
mkdir 1
module load apptainer/1.4.1-gcc-13.3.0-3coysxn

apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven /data/users/theaven/busco_6.1.0--pyhdfd78af_1 busco -i genome.fa -l $Database -m geno -c"$cpu" -f --tar -o 1

cp 1/run*/short_summary.txt ${OutDir}/${OutFile}_short_summary.txt
echo DONE
rm -r $WorkDir
