```bash
du --max-depth=1 --total -h .
```
```bash
salloc --mem=200G
```
```bash
sacct -j 58864242 --format=JobID,JobName,ReqMem,MaxRSS,TotalCPU,AllocCPUS,Elapsed,State,ExitCode,Ntasks,NCPUS,User
sacct -j 19260713 --format=JobID,JobName,ReqMem,MaxRSS,TotalCPU,AllocCPUS,Elapsed,State,ExitCode,Ntasks,NCPUS,User
sacct -j 58565446 --format=JobID,JobName,ReqMem,MaxRSS,TotalCPU,AllocCPUS,Elapsed,State,ExitCode,Ntasks,NCPUS,User
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
find . -type f -name "*.bam"
```
```bash
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 - <<EOF 

EOF

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/earlgrey4.0.6.sif conda activate earlGrey - <<EOF 

EOF
```
```bash
SOURCE_DB_MANIFEST="https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.manifest"
LOCAL_DB="/home/theaven/scratch/apps/fcs/test"
python3 fcs.py db get --mft "$SOURCE_DB_MANIFEST" --dir "$LOCAL_DB/test-only" 
python3 fcs.py db check --mft "$SOURCE_DB_MANIFEST" --dir "$LOCAL_DB/test-only"
srun -p short  --mem 10G --pty bash

SOURCE_DB_MANIFEST="https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest"
LOCAL_DB="/home/theaven/scratch/apps/fcs/all"
python3 fcs.py db get --mft "$SOURCE_DB_MANIFEST" --dir "$LOCAL_DB/gxdb" 

python3 fcs.py db check --mft "$SOURCE_DB_MANIFEST" --dir "$LOCAL_DB/gxdb"

python3 /home/theaven/scratch/apps/fcs/fcs.py screen genome --fasta ./fcsgx_test.fa.gz --out-dir ./gx_out/ --gx-db "/home/theaven/scratch/apps/fcs/all/gxdb/all"  --tax-id 6973

#########################################################################################################################
fasta=~/fcsgx_test.fa.gz
OutDir=~/gx_out
OutFile=testsss
cpu=48

GX_NUM_CORES=${cpu}
export FCS_DEFAULT_IMAGE=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/fcs-gx.sif
python3 ./fcs.py screen genome --fasta ${fasta} --out-dir ${OutDir} --out-basename ${OutFile} --gx-db "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/fcs/gxdb/gxdb" --tax-id 6973 

python3 ./fcs.py db check --mft "$SOURCE_DB_MANIFEST" --dir "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/fcs/gxdb/gxdb"
python3 ./fcs.py clean genome --help

srun -p jic-largemem -c 48 --mem 500G --pty bash
source package 1413a4f0-44e3-4b9d-b6c6-0f5c0048df88
GX_NUM_CORES=48
export FCS_DEFAULT_IMAGE=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/fcs-gx.sif
python3 ./fcs.py screen genome --fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered.fa --out-dir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/fcs/ --gx-db "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/fcs/gxdb/all"  --tax-id 872318




```
```bash
for job in $(squeue -u did23faz | egrep '(AssocMaxJobsLimit|(PartitionTimeLimit))' | awk '{print $1}'); do
    scancel $job
done

for job in $(squeue -u did23faz | egrep 'trim_gal' | awk '{print $1}'); do
    scancel $job
done

for job in $(squeue -u did23faz | grep 'mafft' | awk '{print $1}'); do
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
for job in $(cat logs/pamllog_homN.txt | grep 'Submitted'| awk '{print $4}'); do
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
Delete slurm files that are over 1 day old:
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
singularity overlay create --size 768 ~/overlay.img
sudo singularity shell --overlay ~/overlay.img ~/methylkit1.28.0.sif

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


BiocManager::install("HiCExperiment") #
BiocManager::install("HiContacts")
BiocManager::install("HiCool")
BiocManager::install("HiContactsData")
BiocManager::install("DNAZooData")
BiocManager::install("fourDNData")

install.packages("InteractionSet")
install.packages("ggplot2")
install.packages("GenomicRanges")
install.packages("rtracklayer")
install.packages("dplyr")
install.packages("fourDNData")
install.packages("OHCA")
install.packages("cowplot")
install.packages("purrr")
install.packages("BiocParallel")
install.packages("plyinteractions")
install.packages("multiHiCcompare")
install.packages("tidyr")
install.packages("diffHic")
install.packages("TopDom")
install.packages("GOTHiC")
install.packages("DNAZooData")

BiocManager::install("InteractionSet")
BiocManager::install("ggplot2")
BiocManager::install("GenomicRanges")
BiocManager::install("rtracklayer")
BiocManager::install("dplyr")
BiocManager::install("fourDNData")
BiocManager::install("OHCA")
BiocManager::install("cowplot")
BiocManager::install("purrr")
BiocManager::install("BiocParallel")
BiocManager::install("plyinteractions")
BiocManager::install("multiHiCcompare")
BiocManager::install("tidyr")
BiocManager::install("diffHic")
BiocManager::install("TopDom")
BiocManager::install("GOTHiC")
BiocManager::install("DNAZooData")

library(InteractionSet)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(fourDNData)
library(OHCA)
library(cowplot)
library(purrr)
library(BiocParallel)
library(plyinteractions)
library(multiHiCcompare)
library(tidyr)
library(diffHic)
library(TopDom)
library(GOTHiC)
library(DNAZooData)
```