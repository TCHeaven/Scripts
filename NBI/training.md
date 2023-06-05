## Training perfromed following: http://training.scicomp.jic.ac.uk/docs/cluster_computing_course_book/index.html 
An example of cluster use
```bash
interactive
hostname

mkdir cluster_training
cd cluster_training
#uniprot_sprot.fasta.gz added to folder from https://bit.ly/3BCnzWa
gunzip uniprot_sprot.fasta.gz
source blast+-2.7.1

nano makeblastdb.slurm.sh
	#!/bin/bash
	#SBATCH -p jic-training
	#SBATCH -t 0-00:10

	source blast+-2.7.1
	makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprot

sbatch makeblastdb.slurm.sh #Submitted batch job 54087192

nano test_query.fas
	>query1
	AFQYRPDLLAPYLALPQGEKVQAEYVWIDGDGGLRSKTTTVPKK
	>query2
	YLALPQGEKVQAEYVWIDGDGGLRSKTTTVPKKVTDIGQLRIWDFDGS

nano blast-query.slurm.sh
	#!/bin/bash
	#SBATCH -p jic-training
	#SBATCH -t 0-00:10

	source blast+-2.7.1

	blastp -query test_query.fas -db uniprot -out test.result

sbatch blast-query.slurm.sh #Submitted batch job 54087196
```
Finding software on the HPC
```bash
catalogue
catalogue --production
catalogue --production --search bowtie2
source package bowtie2-2.1.0

catalogue --production --search hisat2
source package hisat2-2.1.0
```
Creating and submitting slurm jobs
```bash
nano build_hisat2_indices.slurm.sh
	#!/bin/bash
	#SBATCH -p jic-short
	#SBATCH --mem 1GB
	#SBATCH -t 0-00:01

	source package hisat2-2.1.0

	srun hisat2-build U00096.3.fasta U00096.3

catalogue --search wget	
source package /tgac/software/testing/bin/wget-1.14
wget https://www.ebi.ac.uk/ena/browser/api/fasta/U00096.3?download=true
mv U00096.3?download=true U00096.3.fasta
sbatch build_hisat2_indices.slurm.sh #sbatch: error: Batch job submission failed: User's group not permitted to use this partition
#with #SBATCH -p jic-training #Submitted batch job 54089330
```
```bash
sinfo -s
#PARTITION           AVAIL  TIMELIMIT   NODES(A/I/O/T)  NODELIST
#jic-training           up    2:00:00         9/3/0/12  n128n[4-9,20-21],n256n[1-4]
#RG-Saskia-Hogenhout    up   infinite          1/0/0/1  j1536n2

# - <institute>-short - performign everyday tasks such as unzipping files
# - <institute>-medium - running quick analyses such as aliging short sequence reads
# - <institute>-long - running long analysis like big phylogenitic bootstrapping jobs

sacct
sacct -j 17346847
sacct -j 17346847 --format=JobID,ReqMem,MaxRSS,TotalCPU #Gives information about how much memory we requested (ReqMem), but also how much we used, MaxRSS as well as the total CPU time used by the job TotalCPU.

whoami
sacct -u did23faz
sacct -j 54089330
sacct -j 54089330 --format=JobID,JobName,MaxRSS,Elapsed

squeue
squeue -u did23faz
squeue -j 54087183

snoderes #A = active, F = free, T = total

nano useless.slurm.sh
	#!/bin/bash
	#SBATCH -p jic-training
	#SBATCH --mem=1G

	sleep 1h

sbatch useless.slurm.sh #Submitted batch job 54089796
squeue -u did23faz
scancel 54089796
```
