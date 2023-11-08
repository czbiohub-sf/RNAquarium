#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Please specify the working directory"
    exit 1
fi

WORKING_DIR="$(realpath "${1}")"

#DEFAULT_ACCESSIONS_LIST="${WORKING_DIR}/data/SRA_accession_list.1.27.23.txt"
#DEFAULT_TOOLS="/hpc/projects/theory_ds/internship/jacob.paras/tools"

DEFAULT_ACCESSIONS_LIST="${WORKING_DIR}/SRA_accession_list.test.txt"
DEFAULT_TOOLS="/hpc/projects/balla_group/sra_experiments/tools"
DEFAULT_ENVNAME="sra"
REMOVEFQ="no"

# click enter on each of these to use the default paths
read -rep $"Enter the path to your accessions list files (press enter to use default: ${DEFAULT_ACCESSIONS_LIST}):"$'\n' accessions_list
ACCESSIONS_LIST=${accessions_list:-${DEFAULT_ACCESSIONS_LIST}}
read -rep $"Enter the directory that has your executables (press enter to use default: ${DEFAULT_TOOLS}):"$'\n' tools
TOOLS=${tools:-${DEFAULT_TOOLS}}
read -rep $"Enter the name of the conda environment to use (press enter to use default: ${DEFAULT_ENVNAME}):"$'\n' envname
ENVNAME=${envname:-${DEFAULT_ENVNAME}}
read -rep $"Remove fastq after STAR run: yes/no (press enter to use default: ${REMOVEFQ}):"$'\n' rmfq
REMOVEFQ=${rmfq:-${REMOVEFQ}}

SLURM_OUTDIR="slurm.out"
echo "removing folder $SLURM_OUTDIR"
rm -rf $SLURM_OUTDIR
mkdir $SLURM_OUTDIR

rm -rf ${WORKING_DIR}/logs/dedup.skip.log
rm -rf ${WORKING_DIR}/logs/dedup.error.log
rm -rf ${WORKING_DIR}/logs/dedup.lengths.log
rm -rf ${WORKING_DIR}/logs/dedup.warning.log

if [ ! -d "${WORKING_DIR}/Dedup_out" ]; then
    mkdir "${WORKING_DIR}/Dedup_out"
fi


sbatch step.5.dedup.sh "${WORKING_DIR}" "${ACCESSIONS_LIST}" "${TOOLS}" "${ENVNAME}" "${REMOVEFQ}"
