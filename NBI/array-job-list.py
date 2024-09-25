#run as python3 array-job-list.py filelist.txt (source python-3.7.2 first)
#script will replace each itteration of the word 'replace' in the template script with each item in the filelist.txt in turn and submit as a separate cluster job
#replace alphafold-template.sh  with name of slurm script to run


import sys, statistics, random, itertools, csv, os, subprocess, time

#For each species in peak file...
file = sys.argv[1]
spec = open(f"{file}",'r').read()
spec= spec.split("\n")

files = []
for i in spec:
    files.append(i)

files = list(set(files))

for i in files:
    if i != "":
        print(i)
        subprocess.Popen([f"cp alphafold-template.sh  {i}.sh; sed -i 's/replace/{i}/g' {i}.sh; sbatch {i}.sh"], shell=True)
        
