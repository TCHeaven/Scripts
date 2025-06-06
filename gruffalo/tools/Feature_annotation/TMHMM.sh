#!/usr/bin/env bash
#SBATCH -J tmhmm
#SBATCH --partition=short
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=8

USAGE='TMHMM.sh <predicted_proteins.fa>'
INFILE=$1
OutDir=$2
echo $USAGE

PRED_PROTEINS=$(echo $INFILE | rev | cut -d "/" -f1 | rev)
CUR_PATH=$PWD
WORK_DIR=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

mkdir -p $WORK_DIR
cd $WORK_DIR

cat $CUR_PATH/$INFILE | tmhmm > tmhmm_out.txt
cat tmhmm_out.txt | grep -v -w 'PredHel=0' > TM_genes_pos.txt
cat tmhmm_out.txt | grep -w 'PredHel=0' > TM_genes_neg.txt

mkdir -p $CUR_PATH/$OutDir
cp -r $WORK_DIR/* $CUR_PATH/$OutDir/.

rm -r $WORK_DIR
