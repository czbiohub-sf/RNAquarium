#!/bin/bash

# Workflow to execute Serratus Lite on a cluster
# A slurm job array is used to search SRA samples in parallel
# The search is performed in 3 steps

# Step 1 make bowtie index
srun --pty --mem 32G bash -l 
module load anaconda
conda activate bowtie2
query=zpoxv.fasta # <---- adjust this when switching query!!!!
VIRALSEQ="/hpc/projects/balla_group/sra_experiments/search_allZfRNAseq/query_virus_data/${query}"
bowtie2-build ${VIRALSEQ} ${VIRALSEQ} 

# Step 2 search
bash start_search.sh

# Step 3 Summarize
conda activate sra
python parse_summarized_files.py --dir summarized