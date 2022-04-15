#!/bin/bash

##### ========================= DESCRIPTION ============================= #####
# A script to index azoreductase gene sequences mined from UHGG for alignment
# 
# Author: Domenick J. Braccia
# Last updated: 12 November 2021
##### =================================================================== #####

## receiving standard input
gene_seqs=$1
gene_index=$2

## run kma index
./external/kma/kma index -i $gene_seqs -o $gene_index
