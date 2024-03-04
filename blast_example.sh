#!/bin/bash

VIRDB=/Data/NCBI_Virus_DB/NCBI_VIRUS

blastn -query $1 -db $VIRDB -outfmt "6 qseqid sseqid evalue pident sskingdom ssciname scomname sblastname stitle" -max_target_seqs 5 -out $1.blast.out -num_threads 20 -perc_identity 30


# After BLAST runs, parse the output

# take top blast hit, use it to find the original contig sequence
cat $1.blast.out | cut -f 1 > positive_ids
grep -A1 -Fwf positive_ids $1 | seqtk seq -A > matches_in_nanopore.fasta


### rename fasta headers using the top blast hits + coverage stats ###

# clean up original nanopore fasta headers
cat matches_in_nanopore.fasta | seqtk seq -A | cut -d " " -f1 > tidy.fasta


# make new headers from blast output
paste <(cat $1.blast.out | cut -f1) <(cat $1.blast.out | cut -f 4) <(cat $1.blast.out | cut -f 9) -d "|" > new_headers

# use while loop to rename fasta headers
while read line; do id=$(echo $line | cut -f1 -d "|"); sed -i "s/$id/$line/" tidy.fasta;done < new_headers


cp tidy.fasta $1.tidy

# add coverage counts to headers
#paste new_headers <(cat $CONTIGS | grep "^>"|cut -d "|" -f4 | uniq -c | sed 's/^[ \t]*//;s/[ \t]*$//' | cut -d " " -f1) > final_headers
#while read line; do id=$(echo $line | cut -f1,2,3,4 -d "|"); sed -i "s/$id/$line/" tidy.fasta;done < final_headers

# get coverage counts (number of times a virus taxon was hit w/ BLAST)
cat new_headers | cut -d "|" -f3 | cut -d " " -f1 | sort | uniq -c > coverage_counts


cat coverage_counts | awk '{ print $1 }' > cov_numbers
cat coverage_counts | awk '{ print $2 }' > cov_ids

FINAL=$1.final.fasta
cat $1 > $FINAL
# add coverage counts to final fasta
while read line;do ID=$(echo $line | cut -d "|" -f 2); sed -i "s/$ID/$line/" $FINAL;done < <(paste -d "|" cov_numbers cov_ids)
