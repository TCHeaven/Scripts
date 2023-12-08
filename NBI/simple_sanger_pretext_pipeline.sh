#!/usr/bin/env bash

#checks input arguments and gives run usage

myInput=''
myDir=''
myHicType=''
myMapDisp=''
mySpecies=''
myDirPath=''
myMapQual=''
myMem=92000
mySort=''
threads=24
myBwa=`which bwa-mem2`
queue="normal"
myHighRes=''

#software identification
mySamtools=`which samtools`
myFilter=''

###checks user
myUser=`whoami`

usage() { 
cat << EOF
Usage: $0 [REQUIRED -i input.fasta, -s species -k kit_type] [OPTIONAL -l,-r, -m mem_request, -q, -d outdir, -t threads, -f,-a queue]

REQUIRED:
	-i [file] Fasta file
	-k [dir]  HiC Kit. Full path to directory containing HiC reads
	-s [text] Species identifier, VGP/TOL format e.g. bGalGal1
OPTIONAL:
	-a [text] name of the LSF queue to use (default $queue)
	-d [file] Full path of user specified output directory
	-m [numeric] Request memory for read mapping (Default $myMem)
	-t [numeric] number of cores to reserve in LSF (default $threads)
	-b	enable legacy bwa mem mode (default is bwa-mem2)
	-f	filter chimeric 5primes (only works with bwa)
	-h	show Usage
	-l	create HiC map with size sorted scaffolds (Default no sort)
	-q	create HiC map with read map quality 10 (Default 0)
	-g	create HiC map in high resolution (Default off)
	-r	soRt the fastafile by id (mutually exlusive with -l)
EOF
	exit 1
}

while getopts "i:d:s:k:m:rblqgfht:a:" opt
   do
     case $opt in
        i ) myInput=$OPTARG
	if [ ! -f $myInput ] ; then echo "File $myInput does not exist"; usage; fi
	;;
        d ) myDirPath=$OPTARG
	if [ ! -d $myDirPath ] ; then echo "Directory $myDirPath does not exist"; usage; fi
	;;
        s ) mySpecies=$OPTARG;;
        k ) myHicType=$OPTARG
	if [ ! -d $myHicType ] ; then echo "Directory $myHicType does not exist"; usage; fi
	;;
    	m ) myMem=$OPTARG;;
	b ) myBwa=`which bwa`;;
	l ) myMapDisp='length';;
    	r ) mySort='1'; myMapDisp='';;
	t ) threads=$OPTARG;;
	a ) queue=$OPTARG;;
	g ) myHighRes='--highRes';;
	q ) myMapQual='10';;
	f ) myFilter=`which filter_five_end`;;
	h | *) usage;;
     esac
done

if [[ ${#myInput} -eq 0 || ${#myHicType} -eq 0 || ${#mySpecies} -eq 0 ]]; then
	usage
fi	 

#check which host we are on; determines bsub variable
echo " "
echo "Checking which host you are running on..."
myHost=`hostname | cut -d "-" -f 1`
myBsubVar="-G grit-grp -q $queue -K"

myFasta=${myInput##*/}
myName=${myFasta%.*}

if [[ ! -f $myInput ]]; then
	echo "$myInput doesn't exist"
	exit 1
fi

if [[ ! -f $myBwa ]];then
	echo "the bwa/bwa-mem2 binary can't be found"
	exit 1
fi

myType=`echo $mySpecies | cut -c 1`

#species type
declare -A speciesType

speciesType=( ["a"]="amphibians" ["b"]="birds" ["c"]="non-vascular-plants" ["d"]="dicots" ["e"]="echinoderms" ["f"]="fish" ["g"]="fungi" ["h"]="platyhelminths" ["i"]="insects" ["j"]="jellyfish" ["k"]="other_chordates" ["l"]="monocotyledons" ["m"]="mammals" ["n"]="nematodes" ["o"]="sponges" ["p"]="protists" ["q"]="arthropods" ["r"]="reptiles" ["s"]="sharks" ["t"]="other_animal_phyla" ["u"]="algae" ["v"]="other_vascular_plants" ["w"]="annelids" ["x"]="molluscs" ["y"]="bacteria" ["z"]="archea" ["-"]="viruses")

mySpeciesType=${speciesType[$myType]}

if [ ${#myDirPath} -eq 0 ]; then
	
	myDirCheck="/lustre/scratch123/tol/teams/grit/$myUser/HiC_Map/"
	
	if [ ! -d $myDirCheck ]; then
			mkdir -p /lustre/scratch123/tol/teams/grit/$myUser/HiC_Map/
			echo "Using default output location"
	else
		echo "Using default output location"
	fi

	myDir="/lustre/scratch123/tol/teams/grit/$myUser/HiC_Map/$mySpeciesType/$myName.$(date '+%Y%m%d:%H:%M:%S')/"

	echo "Creating output directories..."
	mkdir -p $myDir	
	if [ $? -ne 0 ]; then 
		echo "Could not create $myDir"
		exit 1
	fi
else
	myDir="$myDirPath/$myName.$(date '+%Y%m%d:%H:%M:%S')/"
	mkdir -p $myDir
	echo "Using user directory"
fi

echo "Copying fasta..."
cp $myInput $myDir
cd $myDir

echo "Checking fasta header format for bad things....."
sed -i -e '/^>/ s/ .*//' -e  '/^>/ s/|/_/g' -e  '/^>/ s/\./_/g' $myFasta

# in case it needs sorting
if [ "$mySort" == "1" ]; then
        sort_fasta.rb $myFasta > b
	mv b $myFasta
fi

###create bwa index
echo "Creating BWA index..."
bsub $myBsubVar -e index.err -o index.out -M $myMem -R'select[mem>'$myMem'] rusage[mem='$myMem'] span[hosts=1]' "$myBwa index $myFasta"

##map reads
count=1
echo "Mapping Illumina HiC reads..."

mapCpus=$(( threads / 5 ))
CRAM_FILTER=/software/npg/current/bin/cram_filter
for i in $(ls $myHicType/*.cram); do
    CRINDEX=$i.crai
    CFROM=0
    CTO=10000
    NCONTAINERS=`zcat $CRINDEX |wc -l`

    if [ "$NCONTAINERS" == "0" ]; then
	echo "can't find HiC reads in $CRINDEX"
	exit 1
    fi

    while [ "$CFROM" -le "$NCONTAINERS" ]; do
        if [ "$CTO" -gt "$NCONTAINERS" ]; then
          CTO=$NCONTAINERS
        fi
        cmd="$CRAM_FILTER -n $CFROM-$CTO $i - | $mySamtools fastq -@ $mapCpus -F0xB00 -nt - | $myBwa mem -T 10 -t $mapCpus -5SPCp $myFasta - | $mySamtools fixmate -mpu -@ $mapCpus - - | $mySamtools sort -@ $mapCpus --write-index -l9 -o $mySpecies.$count.align.bam -"
        bsub $myBsubVar -e aln_$count.err -o aln_$count.out -n $threads -M $myMem -R'select[mem>'$myMem'] rusage[mem='$myMem'] span[hosts=1]' "$cmd" &		
        echo "Submitting mapping job $count - for $CFROM-$CTO $i"
        ((count++))
        CFROM=$((CTO+1 ))
        CTO=$((CFROM+9999))
    done
done
wait

### check if the output bam files exist, if not die ###
if [[ ! -f "$mySpecies.1.align.bam" ]]; then
	    echo "$mySpecies.1.align.bam doesn't exist. Please check aln_1.err and aln_1.out for the cause and a.e increase memory with -m XYZ"
	    exit 1
fi

###filter BAMs
if [[ -f "$myFilter" ]]; then
	filterCpus=$(( threads / 3 ))
	count=1
	echo "Filtering BAM files..."
	for i in `ls -1 *.align.bam`
		do
			cmd="$mySamtools view -@ $filterCpus -h $i | $myFilter | $mySamtools sort -@ $filterCpus - > $mySpecies.$count.filtered.bam"
			bsub $myBsubVar -e filter_$count.err -o filter_$count.out -n $threads -M $myMem -R'select[mem>'$myMem'] rusage[mem='$myMem'] span[hosts=1]'  "$cmd" &
		echo "Submitting filter job $count"
		((count++))
		done
	wait

### check if the output bam files exist, if not die ###

	if [[ ! -f "$mySpecies.1.filtered.bam" ]]; then
		echo "$mySpecies.1.filtered.bam doesn't exist. Please check filter_1.err and filter_1.out for the cause and a.e increase memory with -m XYZ"
	exit 1
	fi

	rm *.align.bam
fi

if [ `ls -1 *.bam | wc -l` -gt 1 ]; then
	mergeCpus=$(( threads / 2 ))
	echo "Merging BAM files..."
	ls -1 *.bam > bamlist
	cmd="$mySamtools merge -@ $mergeCpus -cpu -b bamlist - |$mySamtools markdup -@ $mergeCpus --write-index - $myName.mapped.bam"
	bsub $myBsubVar -e merge.err -o merge.out -n $threads -M $myMem -R'select[mem>'$myMem'] rusage[mem='$myMem'] span[hosts=1]' "$cmd"
else
	samtools markdup --write-index $mySpecies.1.*.bam $myName.mapped.bam
fi

### check if the output bam files exist, if not die ###
if [[ ! -f "$myName.mapped.bam" ]]; then
	    echo "$myName.mapped.bam doesn't exist. Please check merge.err and merge.out for the cause and a.e increase memory with -m XYZ"
	    exit 1
fi

###create pretext view file
if [ ${#myMapDisp} -eq 0 ]; then 
	myMapDisp="nosort"
fi

if [ ${#myMapQual} -eq 0 ]; then 
	myMapQual="0"
fi

# the hardcoding of the conda binaries is not perfect
pretextMap=/software/grit/conda/envs/curation_v2/bin/PretextMap
bsub $myBsubVar -e pretext.err -o pretext.out -n 5 -M 20000 -R'select[mem>20000] rusage[mem=20000] span[hosts=1]' "$mySamtools view -@ 4 -h $myName.mapped.bam | $pretextMap --sortby $myMapDisp --mapq $myMapQual -o $myName.pretext $myHighRes"

# build a image from pretext
/software/grit/conda/envs/curation_v2/bin/PretextSnapshot -m $myName.pretext --sequences "=full" -c 26 -r 3840 -o .

# cleanup
mkdir logs
mv *.err logs/
mv *.out logs/

rm bamlist *.bam $myFasta* *.csi
rm *mapped.bam

echo "Completed, please check output files in $myDir"

exit 0

