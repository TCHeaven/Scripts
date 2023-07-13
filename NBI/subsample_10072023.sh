#!/bin/bash
#SBATCH --job-name=subsample
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --mem 50G
#SBATCH -c 32
#SBATCH -p jic-long,nbi-long
#SBATCH --time=2-00:00:00

source package /tgac/software/testing/bin/seqkit-0.10.0
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/*/* |grep -v 'S98\|S121\|US1L\|NIC_410R'); do
sample=$(echo $ReadDir | rev | cut -d '/' -f1 | rev)
coverage=$(grep 'mean coverageData' $ReadDir/qualimap/*genome_results.txt | rev | cut -d ' ' -f1 | rev | sed 's@X@@g')
if [[ $coverage < 15 ]]; then
    echo $sample has low coverage: $coverage >> Reports/Low_coverage_report.txt
fi
target_coverage=15
subsampling_ratio=$(bc <<< "scale=2; $target_coverage / $coverage")
for reads in $(ls $ReadDir/*.fq.gz); do
OutFile=$(basename $reads | cut -d '.' -f1)_subsampled.fq.gz
OutDir=${ReadDir}/subsampled
OutDir2=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/M_persicae/${sample}
mkdir $OutDir
mkdir $OutDir2
echo $sample raw reads have average coverage of ${coverage}, these were subsampled with ratio $subsampling_ratio >> Reports/Subsampling_report.txt
seqkit sample -s 1234 -p $subsampling_ratio -o ${OutDir2}/${OutFile} $reads
ln -s ${OutDir2}/${OutFile} ${OutDir}/${OutFile}
done
if [ -e ${OutDir2}/${OutFile} ]; then
    echo ${sample} Done
rm -r /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/M_persicae_Bass_SRA/$sample
else 
rm ${OutDir}/*_subsampled.fq.gz
for reads in $(ls $ReadDir/*.fq.gz); do
name=$(echo $reads | rev | cut -d '/' -f1 | rev)
cp $reads ${OutDir2}/${name}
ln -s ${OutDir2}/${name} ${OutDir}/${name}
rm -r /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/M_persicae_Bass_SRA/$sample
done
fi
done