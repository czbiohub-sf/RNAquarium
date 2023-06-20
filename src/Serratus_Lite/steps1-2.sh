#!/bin/bash

# Workflow to execute Serratus Lite on a cluster
# A slurm job array is used to search SRA samples in parallel
# The search is performed in 3 steps

if [ $# -ne 1 ]; then
    echo "Please specify the parent directory of the query virus data"
    exit 1
fi

# this script should be a sibling file to the directory that contains the virus sequences
basedir="$(realpath "${1}")"

DEFAULT_ACCESSIONS_LIST="/hpc/scratch/group.theory/jparas/Zebrafish-RNA-Quarium/src/nonhost/data/SRA_accession_list.1.27.23.txt"
DEFAULT_SEQBASEDIR="/hpc/projects/theory_ds/internship/jacob.paras/Gsnap_out"
DEFAULT_TOOLS="/hpc/projects/theory_ds/internship/jacob.paras/tools"

# click enter on each of these to use the default paths
read -rep $"Enter the path to your accessions list files (press enter to use default):"$'\n' accessions_list
accessions_list=${accessions_list:-${DEFAULT_ACCESSIONS_LIST}}
read -rep $"Enter the path to your RNASeq output directory (press enter to use default):"$'\n' seqbasedir
seqbasedir=${seqbasedir:-${DEFAULT_SEQBASEDIR}}
read -rep $"Enter the directory that has your executables (press enter to use default):"$'\n' tools
tools=${tools:-${DEFAULT_TOOLS}}

# edit this path to point to your installation directory
bowtie2="${tools}/bowtie2-2.4.5-linux-x86_64"

if [ ! -d "${basedir}" ]; then
    echo "Directory ${basedir} does not exist"
    exit 1
fi

# Step 1 make bowtie index
for virus in "${basedir}"/query_virus_data/*.fa ; do
    echo "Running query on ${virus}..."
    "${bowtie2}/bowtie2-build" "${virus}" "${virus}"
    virus_name="$(basename "${virus%.*}")"
    res_folder=${basedir}/results/${virus_name}
    if [ -e "$res_folder" ]; then
        rm -rf "$res_folder"
    fi    
    mkdir -p "$res_folder"
    mkdir "$res_folder/summarized"
    mkdir "$res_folder/bam"
done

# Step 2 run search on hpc using job array
if [ -e slurm.out/ ] ; then rm -rf slurm.out/ ; fi
mkdir -p slurm.out
cd "${basedir}/scripts"
sbatch hpc_search.sh "${basedir}" "$(realpath "${accessions_list}")" "$(realpath "${seqbasedir}")" "$(realpath "${tools}")"
