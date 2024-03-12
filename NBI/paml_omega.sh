#!/bin/bash

mapfile -t Seqfiles < Seqfiles.txt
mapfile -t TreeFiles < Treefiles.txt
mapfile -t OutDirs < Outdirs.txt
mapfile -t OutFiles < Outfiles.txt

Seqfile=${Seqfiles[${SLURM_PROCID}]}
TreeFile=${TreeFiles[${SLURM_PROCID}]}
OutDir=${OutDirs[${SLURM_PROCID}]}
OutFile=${OutFiles[${SLURM_PROCID}]}

mkdir -p ${OutDir}
mkdir workdir_${SLURM_PROCID}
cd workdir_${SLURM_PROCID}

echo seqfile = $Seqfile > ${SLURM_PROCID}_codeml.ctl
echo treefile = $TreeFile >> ${SLURM_PROCID}_codeml.ctl
echo outfile = ${SLURM_PROCID}_results.out >> ${SLURM_PROCID}_codeml.ctl
echo noisy = 0 >> ${SLURM_PROCID}_codeml.ctl
echo verbose = 0 >> ${SLURM_PROCID}_codeml.ctl
echo runmode = 0 >> ${SLURM_PROCID}_codeml.ctl
echo seqtype = 1 >> ${SLURM_PROCID}_codeml.ctl
#echo cleandata = 1 >> ${SLURM_PROCID}_codeml.ctl #This will remove sites with ambiguity data, eg. 'N', default is =0
#echo CodonFreq = 7 >> ${SLURM_PROCID}_codeml.ctl #This is more accurate but will require more resources/time
echo CodonFreq = 2 >> ${SLURM_PROCID}_codeml.ctl
echo clock = 0 >> ${SLURM_PROCID}_codeml.ctl
echo model = 0 >> ${SLURM_PROCID}_codeml.ctl
echo NSsites = 0 >> ${SLURM_PROCID}_codeml.ctl
echo icode = 0 >> ${SLURM_PROCID}_codeml.ctl
echo fix_kappa = 0 >> ${SLURM_PROCID}_codeml.ctl
echo kappa = 2 >> ${SLURM_PROCID}_codeml.ctl
echo fix_omega = 0 >> ${SLURM_PROCID}_codeml.ctl
echo omega = 0.4 >> ${SLURM_PROCID}_codeml.ctl
echo ncatG = 8 >> ${SLURM_PROCID}_codeml.ctl

codeml ${SLURM_PROCID}_codeml.ctl 2>&1 >> ${SLURM_PROCID}_stop.txt

cp ${SLURM_PROCID}_results.out ${OutDir}/${OutFile}
cp ${SLURM_PROCID}_stop.txt ${OutDir}/${OutFile}_stop.txt

cd ..