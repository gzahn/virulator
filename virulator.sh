#!/bin/bash

set -ueo pipefail

# To run: plant_virus_blast "host sci_name" /path/to/Nanopore/fastqs

# Depends: blastn, minimap2, samtools, flye 2.9+, R 4+, seqtk,
# external files at:
# Zahn, Geoffrey. (2022). NCBI Virus BLAST Database [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7250322


echo "Usage: bash virulator.sh 'Host sciname' path/to/nanopore/reads/directory/ '[assembly|raw]'"

# check for host name argument
if [ -z "$1" ]
  then
    echo "No host name specified."
    exit 1
fi

# check for path argument
if [ -z "$2" ]
  then
    echo "No path/to/nanopore/reads/directory/ was specified."
    exit 1
fi

# check for "assembly" argument
if [ -z "$3" ]
  then
    echo "Need to specify either assembly or raw."
    echo " "
    echo "e.g.,: bash virulator.sh 'brassica oleracea' ./broccoli 'assembly'"
    exit 1
fi


# Set max cores as avail - 1
NCORES="$(nproc --all | awk '{ncor = $1 - 1; print ncor}')"
echo "$NCORES cores will be used"

# download reference genome (datasets)

# look for datasets software in path
# if not found, then download and make executable

#If downloaded it makes it executable.
if ! command -v datasets &> /dev/null
then
  curl -o datasets "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets"
  chmod +x dataset
  echo "Datasets downloaded"
  export PATH=$PATH:$(pwd)
  echo "datasets temporarily added to PATH"
else
  echo "datasets already in PATH, skipping download."
fi


# check for bbtools
if ! command -v bbmap.sh &> /dev/null
then
  echo "bbtools does not appear to be installed in your PATH, idiot."
exit 1
fi


######## Check for NCBI virus database #########

# check for pre-existing virus database in home dir

FILE=~/NCBI_VIRUS_DB/NCBI_VIRUS.nto
DB_MD5=$(md5sum $FILE | cut -f1 -d " ")
if [ $DB_MD5 == "5f98018930d12c35f042e3030dcb4830" ]; then
    echo "$FILE exists and md5sum is correct...skipping virus database download."
else

# download virus database
wget https://zenodo.org/record/7255748/files/NCBI_VIRUS.ndb
wget https://zenodo.org/record/7255748/files/NCBI_VIRUS.nhr
wget https://zenodo.org/record/7255748/files/NCBI_VIRUS.nin
wget https://zenodo.org/record/7255748/files/NCBI_VIRUS.not
wget https://zenodo.org/record/7255748/files/NCBI_VIRUS.nsq
wget https://zenodo.org/record/7255748/files/NCBI_VIRUS.ntf
wget https://zenodo.org/record/7255748/files/NCBI_VIRUS.nto

# move virus database into home directory
mkdir ~/NCBI_VIRUS_DB
mv NCBI_VIRUS.* ~/NCBI_VIRUS_DB

# re-try hash of file
DB_MD5=$(md5sum $FILE | cut -f1 -d " ")
if [ $DB_MD5 != "5f98018930d12c35f042e3030dcb4830" ]; then
 echo "Your download of the virus database was unsuccessful or corrupted. Exiting..."
 exit 1
fi

fi

# start documenting versions, etc
datasets --version > documentation.txt
echo "Run on:" >> documentation.txt
date >> documentation.txt
blastn -version >> documentation.txt
echo "BBTools version:"
bbversion.sh >> documentation.txt
samtools version | head -2 >> documentation.txt
spades.py --version >> documentation.txt

# Set sirus database variable
VIRDB=$HOME/NCBI_VIRUS_DB/NCBI_VIRUS


# make new name for genome folder
HOST=${1/ /_}
mkdir $HOST


# download the host genome
echo "Downloading specified host genome assembly..."
datasets download genome taxon "$1" --reference --include genome --filename $HOST.zip
echo "$1"
echo "Host genome assembly downloaded."

# unzip host genome
unzip -o $HOST.zip -d $HOST

# pull all genomic fna files into a single file
find $HOST -name "*.fna" -type f | xargs cat > $HOST/full_host_genome.fna

# cleanup redundant fastas
rm -r $HOST/ncbi_dataset



# map reads to reference (bbmap)
REFERENCE=$HOST/full_host_genome.fna
QUERY_DIR=$2
# make concatenated query file
zcat $QUERY_DIR/*.fastq.gz > ./complete_nanopore_run.fastq
QUERY=./complete_nanopore_run.fastq

# map reads and build sam file ###########

# convert to fasta and clean up headers
cat $QUERY | sed 's/ .*//' > clean.fastq

echo "Running mapPacBio.sh to align reads to host genome..."

# run bbmap
mapPacBio.sh ref=$REFERENCE in=clean.fastq out=hostmap.sam nodisk qin=33 covstats=host_coverage.stats

# keep only reads that didn't map well
samtools view -f 4 hostmap.sam | samtools fasta > unmapped_reads.fasta
samtools view -f 4 hostmap.sam | samtools fastq > unmapped_reads.fastq

# if user specifies 'assembly' then run assembly before BLAST #######################
if [ $3 == "assembly" ]; then

 echo "'assembly' specified. Attempting contig assembly before BLAST."

# build contigs (flye / Spades)
# if long enough >3000, use flye

echo "Your N50 for raw reads is:"

cat unmapped_reads.fastq | seqtk seq -A > unmapped_reads.fasta
N50=$(stats.sh -in=unmapped_reads.fasta | grep "Main genome contig N/L50" | cut -d "/" -f3)
echo $N50
# set minimum read N50 for Flye to run
min_flye=3000


if [ "$N50" -gt "$min_flye" ]; then
 # Flye assembly
 echo "Using Flye for assembly"
 flye --nano-raw unmapped_reads.fasta -o Flye_Assembly --meta 
else
 # spades assembly
 echo "Using Spades for assembly"
 spades.py -t 16 -o Spades_Assembly -s unmapped_reads.fastq
fi


# map raw reads onto contigs to get coverage stats for each contig (using bbmap)


# blast against known viruses (blastn) # permissive/sloppy


# make variable for contigs
# if flye was used, then use flye path to contigs
# else, if spades was used used that path
# else, if neither worked....use raw reads

# define possible outputs
SPADES_OUTPUT=./Spades_Assembly/scaffolds.fasta
FLYE_OUTPUT=./Flye_Assembly/assembly.fasta

if test -f "$FLYE_OUTPUT"; then
CONTIGS="$FLYE_OUTPUT";
echo "Using Flye contigs for BLAST";
elif test -f "$SPADES_OUTPUT"; then
CONTIGS="$SPADES_OUTPUT";
echo "Using Spades contigs for BLAST";
else
CONTIGS=unmapped_reads.fasta;
echo "Neither Flye nor Spades worked. Using raw non-host reads for BLAST."
fi



# run blast on contigs (whichever worked, or just unmapped reads)
blastn -query $CONTIGS -db $VIRDB -outfmt "6 qseqid sseqid evalue pident sskingdom ssciname scomname sblastname stitle" -max_target_seqs 1 -out blast.out -num_threads $NCORES


# If user specifies "raw" then skip the whole assembly attempt
elif [ $3 == "raw" ]; then
echo "'raw' specified. Running BLAST on unmapped reads (skipping assembly attempt)"

# convert unmapped_reads.fq to fasta
cat unmapped_reads.fastq | seqtk seq -A > unmapped_reads.fasta

# set "contigs" to be unmapped raw reads
CONTIGS=unmapped_reads.fasta
blastn -query $CONTIGS -db $VIRDB -outfmt "6 qseqid sseqid evalue pident sskingdom ssciname scomname sblastname stitle" -max_target_seqs 1 -out blast.out -num_threads $NCORES

fi ###################################################################################

# give user message on BLAST status, or abort if unsuccessful

if test -f blast.out; then
  echo "BLAST successfully completed.";
else
  echo "BLAST did not complete successfully";
  exit 1
fi

# After BLAST runs, parse the output

# take top blast hit, use it to find the original contig sequence
cat blast.out | cut -f 1 > positive_ids
grep -A1 -Fwf positive_ids $CONTIGS | seqtk seq -A > matches_in_nanopore.fasta


### rename fasta headers using the top blast hits + coverage stats ###

# clean up original nanopore fasta headers
cat matches_in_nanopore.fasta | seqtk seq -A | cut -d " " -f1 > tidy.fasta


# make new headers from blast output
# Still need to add coverage stats to header as well ................
paste <(cat blast.out | cut -f1) <(cat blast.out | cut -f 4) <(cat blast.out | cut -f 9) -d "|" > new_headers

# use while loop to rename fasta headers
while read line; do id=$(echo $line | cut -f1 -d "|"); sed -i "s/$id/$line/" tidy.fasta;done < new_headers

# add coverage counts to headers
#paste new_headers <(cat $CONTIGS | grep "^>"|cut -d "|" -f4 | uniq -c | sed 's/^[ \t]*//;s/[ \t]*$//' | cut -d " " -f1) > final_headers
#while read line; do id=$(echo $line | cut -f1,2,3,4 -d "|"); sed -i "s/$id/$line/" tidy.fasta;done < final_headers

# get coverage counts (number of times a virus taxon was hit w/ BLAST)
cat new_headers | cut -d "|" -f3 | cut -d " " -f1 | sort | uniq -c > coverage_counts
cat coverage_counts | awk '{ print $1 }' > cov_numbers
cat coverage_counts | awk '{ print $2 }' > cov_ids



# rename final output file

if [ $3 == "raw" ]; then
FINAL="$HOST.raw.final.fasta"
mv tidy.fasta $FINAL
elif [ $3 == "assembly" ]; then
FINAL="$HOST.assembled.final.fasta"
mv tidy.fasta $FINAL
fi

# cleanup
if [ $3 == "raw" ]; then
OUTDIR="$HOST.raw.output";
elif [ $3 == "assembly" ]; then
OUTDIR="$HOST.assembled.output";
fi

# move final files to new directory named after host
mkdir $OUTDIR
mv blast.out clean.fastq documentation.txt final_headers host_coverage.stats hostmap.sam matches_in_nanopore.fasta new_headers positive_ids unmapped_reads.fast* $OUTDIR
rm -r $HOST $HOST.zip complete_nanopore_run.fastq
echo "Final output files are in $OUTDIR"
mv $FINAL $OUTDIR

# add coverage counts to final fasta
while read line;do ID=$(echo $line | cut -d "|" -f 2); sed -i "s/$ID/$line/" $OUTDIR/$FINAL;done < <(paste -d "|" cov_numbers cov_ids)


# remove assembly directory, if present
if test -d Spades_Assembly; then
mv $SPADES_OUTPUT $OUTDIR
echo "Spades assembly file is named $SPADES_OUTPUT"
rm -r Spades_Assembly;

elif test -d Flye_Assembly; then
mv $FLYE_OUTPUT $OUTDIR
echo "Flye assembly file is named $FLYE_OUTPUT"
rm -r Flye_Assembly;
fi


echo "Fasta with annotated headers is named $FINAL"
echo "If assembly was attempted and $FINAL is empty, check the assembly file."
