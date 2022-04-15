#!/bin/bash
#SBATCH -J quant_hmp2_array # Job name 
#SBATCH -o logs/quant_hmp2_array_%j_%a.o # output file
#SBATCH -e logs/quant_hmp2_array_%j_%a.e # error file
#SBATCH --mail-user=dbraccia@umd.edu # Email for job info
#SBATCH --mail-type=fail,end # Get email for end, and fail
#SBATCH --time=18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --array=1-708

# MAKE SURE THIS COMMAND GETS RUN IN A TERMINAL BEFORE RUNNING THIS SCRIPT
#ls -1a /fs/cbcb-lab/hall/data/hmp2/mgx-hmp2/*_1.fq.gz >> processed/hmp2_readfile_paths.txt

echo "- loading in modules"
module purge
module load salmon

file=`head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1`
echo "the current file path is: $file"
azored_genes_index=$2
hmp2_gene_quants=$3

current_sample_path=${file/_1.fq.gz/} 
echo "- the current sample path is $current_sample_path"

current_sample_name=`basename $current_sample_path`
echo "- Processing sample $current_sample_name"

# making quants dir if not already
mkdir -p processed/salmon_out/

echo "- Getting quants for: $current_sample_name"
time salmon quant -i $azored_genes_index -l A \
         -1 ${current_sample_path}_1.fq.gz \
         -2 ${current_sample_path}_2.fq.gz \
         -p 16 --validateMappings -o ${hmp2_gene_quants}/${current_sample_name}

echo "- DONE"