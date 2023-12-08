Setting up github on the Norwich biosciences HPC.
```bash
interactive
source package /tgac/software/production/bin/git-latest
mkdir git_repos
cd git_repos

#Each of the following are a separate repository
git clone https://github.com/TCHeaven/Aphididae
git clone https://github.com/TCHeaven/Psyllidae
git clone https://github.com/TCHeaven/Wrappers

#NOTE: cannot push to github from the head node
cd /hpc-home/did23faz/git_repos/Aphididae
git add .
git commit . -m "" 
git push origin main

cd /hpc-home/did23faz/git_repos/Psyllidae
git add .
git commit . -m ""
git push origin main

cd /hpc-home/did23faz/git_repos/Wrappers/NBI
git add .
git commit . -m ""
git push origin main

cd /hpc-home/did23faz/git_repos/Scripts/NBI
git add .
git commit . -m ""
git push origin main
```
```bash
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
```
