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

sbatch step.1.fastqDump.sh "${WORKING_DIR}" "${ACCESSIONS_LIST}" "${TOOLS}"
