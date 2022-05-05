#!/bin/bash

# DESCRIPTION: a script to extract nt sequences of all significant hits for
# 	selected genome (standard input $1). No other user input necessary
#
# simply run: `bash code/cluster/mine_nt_seqs.sh MGYG-HGUT-00236` for example.
#
#
# Author: Domenick J. Braccia
# Date: 5 May 2022

## import modules
module purge
module load bedtools

## import data from standard input
genome_id=$1
sighits=results/from-GutFunFind/Azoreductase.sighits_hmm.tsv

## generate bed file for seqs of interest

parent_fldr=${genome_id::-2}
genome_fldr=../../data/UHGG/uhgg_catalogue/$parent_fldr/$genome_id/genome

## search genome_id in sighits file to get all hits
seq_ids=$(grep "$genome_id" $sighits | awk ' { print $2 } ')
gene_names=$(grep "$genome_id" $sighits | awk ' { print $6 } ')
echo "-- SET OF seq_ids: " $seq_ids
echo "-- SET OF gene_names: " $gene_names

bed_out=./results/mined_seqs/${genome_id}.bed

if [ -e $bed_out ]; then
	rm -r $bed_out
	echo '-- removed old bed_out file'
fi

for id in $seq_ids; do
	grep "$id" $genome_fldr/${genome_id}.gff | \
		awk ' OFS=" " {print $1"\t", $4"\t", $5"\t" } ' | tr -d " " | tr '\n' '\t' >> $bed_out 
	grep "$id" $sighits | awk ' {print $6} ' >> $bed_out
done

## search genome file using bed file generated in previous step

fasta_in=../../data/UHGG/uhgg_catalogue/$parent_fldr/$genome_id/genome/${genome_id}.fna
echo "-- BED out FILE IS: " $bed_out
fasta_out=./results/mined_seqs/${genome_id}.fna
echo "-- FASTA out FILE IS: " $fasta_out
bedtools getfasta -fi $fasta_in -bed $bed_out -fo $fasta_out -split


