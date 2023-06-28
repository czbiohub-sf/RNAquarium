#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Please specify the working directory"
    exit 1
fi

WORKING_DIR="$(realpath "${1}")"

DEFAULT_ACCESSIONS_LIST="${WORKING_DIR}/data/SRA_accession_list.1.27.23.txt"
DEFAULT_TOOLS="/hpc/projects/theory_ds/internship/jacob.paras/tools"

# click enter on each of these to use the default paths
read -rep $"Enter the path to your accessions list files (press enter to use default ${DEFAULT_ACCESSIONS_LIST}):"$'\n' accessions_list
ACCESSIONS_LIST=${accessions_list:-${DEFAULT_ACCESSIONS_LIST}}
read -rep $"Enter the directory that has your executables (press enter to use default ${DEFAULT_TOOLS}):"$'\n' tools
TOOLS=${tools:-${DEFAULT_TOOLS}}

SLURM_OUTDIR="slurm.out"

echo "removing folder $SLURM_OUTDIR"
rm -rf $SLURM_OUTDIR
mkdir $SLURM_OUTDIR

echo "removing folder fastq"
rm -rf fastq
mkdir fastq

echo "removing folder prefetched"
rm -rf prefetched
mkdir prefetched

echo "checking logs directory"
if [ ! -d "${WORKING_DIR}/logs" ]
    then
        mkdir "${WORKING_DIR}/logs"
fi

echo "removing previous log files"

if [ -e $log0 ]
    then
        rm $log0
fi
if [ -e $log1 ]
    then
        rm $log1
fi
if [ -e $log2 ]
    then
        rm $log2
fi

sbatch step.1.fastqDump.sh "${WORKING_DIR}" "${ACCESSIONS_LIST}" "${TOOLS}"
