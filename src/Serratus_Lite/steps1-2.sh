#!/bin/bash

# Workflow to execute Serratus Lite on a cluster
# A slurm job array is used to search SRA samples in parallel
# The search is performed in 3 steps

if [ $# -ne 1 ]; then
    echo "Please specify the parent directory of the query virus data"
fi

# this script should be a sibling file to the directory that contains the virus sequences
basedir="$(realpath "${1}")"
bowtie2="/hpc/projects/theory_ds/internship/jacob.paras/tools/bowtie2-2.4.5-linux-x86_64"

# Step 1 make bowtie index
for virus in $(find "${basedir}/query_virus_data" -maxdepth 1 -type f)
do
    echo "Running query on ${virus}..."
    "${bowtie2}/bowtie2-build" "${virus}" "${virus}"
    virus_name="$(basename "${virus%.*}")"
    res_folder=${basedir}/results/${virus_name}
    if [ -e "$res_folder" ]; then
        rm -rf "$res_folder"
    fi    
    mkdir -p "$res_folder"
    mkdir "$res_folder/slurm.out"
    mkdir "$res_folder/summarized"
    mkdir "$res_folder/bam"
done

# Step 2 run search on hpc using job array
# sbatch search_nonZFv2.sh 
