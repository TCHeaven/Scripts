```bash
sacct -j 54241155 --format=JobID,JobName,ReqMem,MaxRSS,TotalCPU,AllocCPUS,Elapsed,State,ExitCode
```
```bash
Jobs=$(squeue -u did23faz| grep 'omniHiC'  | wc -l)
echo x
while [ $Jobs -gt 5 ]; do
    sleep 300s
    printf "."
    Jobs=$(squeue -u did23faz| grep 'omniHiC'  | wc -l)
done

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/yahs/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
```
```bash
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
```
```bash
for job in $(squeue -u did23faz | grep 'paml' | awk '{print $1}'); do
pattern=$(grep -A 1 'Seqfile:' slurm.${job}.out | sed 's@Seqfile:@@g' | grep -v '^$')
grep -v $pattern temp_files.txt > temp_files2.txt
cat temp_files2.txt > temp_files.txt
done

Jobs=$(squeue -u did23faz| grep 'paml'  | wc -l)
for line in $(cat temp_files.txt); do
echo x
echo $line
while [ $Jobs -gt 123 ]; do
    sleep 300s
    printf "."
    Jobs=$(squeue -u did23faz| grep 'paml'  | wc -l)
done

#loop through files in directory and run ask.sh on each of them
directory=
for 


```
