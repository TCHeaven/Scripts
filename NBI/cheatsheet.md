```bash
du --max-depth=1 --total -h .
```
```bash
salloc --mem=200G
```
```bash
sacct -j 57795761 --format=JobID,JobName,ReqMem,MaxRSS,TotalCPU,AllocCPUS,Elapsed,State,ExitCode
57715532
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
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 - <<EOF 

EOF
```
```bash
for job in $(squeue -u did23faz | egrep '(AssocMaxJobsLimit|(PartitionTimeLimit))' | awk '{print $1}'); do
    scancel $job
done

for job in $(squeue -u did23faz | grep 'bsmap' | awk '{print $1}'); do
    scancel $job
done

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
```
```bash
for job in $(cat logs/pamllog.txt | grep 'Submitted'| awk '{print $4}'); do
sacct -j "$job" --format=JobID,JobName,ReqMem,MaxRSS,TotalCPU,AllocCPUS,Elapsed,State,ExitCode | grep 'paml' |  awk '{print $6}' >> temp_times.txt
done

for job in $(echo {57795755..57795777}); do
x=$(sacct -j "$job" --format=JobID,JobName,ReqMem,MaxRSS,TotalCPU,AllocCPUS,Elapsed,State,ExitCode | grep 'COMPLETED'| grep 'bash'| awk '{print $1}')
grep 'Inchworm file:' slurm.${x}.out
done

Inchworm file:

57758418-65

mapfile -t times < temp_times.txt

time_to_seconds() {
    local h=$(echo $1 | cut -d: -f1 | sed 's/^0*//')  # Remove leading zeros
    local m=$(echo $1 | cut -d: -f2 | sed 's/^0*//')
    local s=$(echo $1 | cut -d: -f3 | sed 's/^0*//')
    echo $((10#h * 3600 + 10#m * 60 + 10#s))
}

total_seconds=0
for time in "${times[@]}"; do
    total_seconds=$((total_seconds + $(time_to_seconds "$time")))
done

average_seconds=$((total_seconds / ${#times[@]}))

average_time=$(printf "%02d:%02d:%02d\n" $((average_seconds / 3600)) $(((average_seconds % 3600) / 60)) $((average_seconds % 60)))

echo "Average time: $average_time"
```
Delete temp directories:
```bash
running_jobs=$(squeue -u did23faz | awk 'NR>1{print $1}')
for folder in tmp_*; do
    # Extract the job ID from the folder name
    folder_job_id=$(echo "$folder" | awk -F "_" '{print $2}')
    
    # Check if the folder's job ID is not in the list of running jobs
    if [[ -n "$folder_job_id" && ! "$running_jobs" =~ "$folder_job_id" ]]; then
        echo "Deleting folder: $folder"
        rm -rf "$folder"
    else
        echo "Keeping folder: $folder"
    fi
done

tmp_57221669
Keeping folder: tmp_57221674
Deleting folder: tmp_57291096

```
```bash
# Define the script to run
cleanup_script="for file in \$(ls dna_qc/Myzus/persicae/wouters/*/subsampled/trim_galore/*.fq.gz dna_qc/Myzus/persicae/singh/*/subsampled/trim_galore/*.fq.gz); do rm \$file; done"

# Run the script every 5 minutes for 10 hours
for ((i=0; i<120; i++)); do
    echo "Running cleanup script (iteration $((i+1)))"
    eval "$cleanup_script"
    sleep 300s  # 300 seconds = 5 minutes
done

ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/singh/*/subsampled/trim_galore/*.fq.gz
```
```bash
find . -name "slurm*" -ctime +0 -print > temp_delete_slurm.txt

for file in $(cat temp_delete_slurm.txt); do
    rm $file
done
```
```bash
tar -xzvf folder_to_decompress.tar.gz
tar -czvf archive_name.tar.gz folder_to_compress
```
```bash
cd ~
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/methylkit1.28.0.sif .
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/genomation1.34.0.sif .
singularity overlay create --size 768 ~/overlay.img
sudo singularity shell --overlay ~/overlays/overlay.img ~/overlays/genomation1.34.0.sif

R --vanilla
setrepo = getOption("repos")
setrepo["CRAN"] = "http://cran.uk.r-project.org"
options(repos = setrepo)
rm(setrepo)
install.packages("lattice")
install.packages("ggplot2")
install.packages("gplots")
Sys.setenv(RGL_USE_NULL = TRUE)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("methylKit")
q()
```