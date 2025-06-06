#!/bin/bash
#SBATCH -J interproscan
#SBATCH --partition=short
#SBATCH --mem=20G
#SBATCH --cpus-per-task=16

# run_interproscan.sh
# The current version of interproscan only works with Java version 11
# Note - the latest version of interproscan doesnt work on the cluster
# because of GLIBC cant be updated.

CUR_PATH=$PWD
IN_FILE=$1
ORGANISM=$2
STRAIN=$3
IN_NAME=$(basename $IN_FILE)

WORK_DIR=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}

mkdir -p $WORK_DIR
cd $WORK_DIR
cp $CUR_PATH/$IN_FILE $IN_NAME
sed -i -r 's/\*/X/g' $IN_NAME

# /mnt/shared/scratch/jconnell/apps/interproscan-5.55-88.0/interproscan.sh -appl CDD,COILS,Gene3D,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,SFLD,SMART,SUPERFAMILY,TIGRFAM -goterms -iprlookup -pa -i $IN_NAME
#/mnt/shared/scratch/theaven/apps/conda/envs/interproscan/share/InterProScan/interproscan.sh -appl CDD,COILS,Gene3D,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,SFLD,SMART,SUPERFAMILY,TIGRFAM -goterms -iprlookup -pa -i $IN_NAME
/home/theaven/scratch/apps/interproscan/interproscan-5.57-90.0/interproscan.sh -appl CDD,COILS,Gene3D,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,SFLD,SMART,SUPERFAMILY,TIGRFAM -goterms -iprlookup -pa -i $IN_NAME


OUT_DIR=$CUR_PATH/$4/raw
mkdir -p $OUT_DIR
cp *.gff3 $OUT_DIR/.
cp *.tsv $OUT_DIR/.
cp *.xml $OUT_DIR/.

rm -r $WORK_DIR
