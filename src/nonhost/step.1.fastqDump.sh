#!/bin/bash

#SBATCH --job-name=fastqDump
#SBATCH --time=14-00:00:00
#SBATCH --array=1-54189%100
#SBATCH --partition preempted
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH -e slurm.out/slurm-%A_%a.err
#SBATCH -o slurm.out/slurm-%A_%a.out

#setting directories
WORKING_DIR="${1}"
ACCESSIONS_LIST="${2}"
TOOLS="${3}"
SRA_BIN="${TOOLS}/sratoolkit.3.0.0-ubuntu64/bin"

declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))

module load anaconda
conda activate zf_pipeline

# declare arrays
readarray -t ACCESSIONS < <(cat "${ACCESSIONS_LIST}") #54189 64G 4cpu

#output directories
FDIR=${WORKING_DIR}/fastq
PDIR=${WORKING_DIR}/prefetched
STDDIR=${WORKING_DIR}/stdout_stderr

echo "accession: ${ACCESSIONS[$idx]}"

############################
# proceed with current job #
############################
echo "processing..." >> ${WORKING_DIR}/logs/STAR.process.log
#sleep 0-120 seconds to space out 100 concurrent jobs from doing prefetch at the same time (prevent timing outs)
sleep $((idx % 120))
cd $WORKING_DIR

###############
# fasterqdump #
###############
#create folder for fastq
if [ -e ${FDIR}/${ACCESSIONS[$idx]} ]
then
    rm -rf ${FDIR}/${ACCESSIONS[$idx]}
fi
mkdir -p ${FDIR}/${ACCESSIONS[$idx]}

if [ -e ${STDDIR}/${ACCESSIONS[$idx]} ]
then
    rm -rf ${STDDIR}/${ACCESSIONS[$idx]}
fi
mkdir -p ${STDDIR}/${ACCESSIONS[$idx]}

#prefetch
${SRA_BIN}/prefetch --max-size 1t --force ALL --output-directory ${PDIR} ${ACCESSIONS[$idx]} 1> ${STDDIR}/${ACCESSIONS[$idx]}/prefetch.stdout.txt 2> ${STDDIR}/${ACCESSIONS[$idx]}/prefetch.stderr.txt

#fasterq-dump
cd ${PDIR}
${SRA_BIN}/fasterq-dump --split-3 --mem 24G --outdir ${FDIR}/${ACCESSIONS[$idx]} ${ACCESSIONS[$idx]} 1> ${STDDIR}/${ACCESSIONS[$idx]}/fqdump.stdout.txt 2> ${STDDIR}/${ACCESSIONS[$idx]}/fqdump.stderr.txt  && {
    echo fasterq-dump: no error
} || {
    echo fasterq-dump encounterd error, reverting back to using fastqdump
    ${SRA_BIN}/fastq-dump --split-3 --disable-multithreading --outdir ${FDIR}/${ACCESSIONS[$idx]} ${ACCESSIONS[$idx]} 1> ${STDDIR}/${ACCESSIONS[$idx]}/fqdump.stdout.txt 2> ${STDDIR}/${ACCESSIONS[$idx]}/fqdump.stderr.txt
}

#remove prefretched data
cd ..
rm -rf ${PDIR}/${ACCESSIONS[$idx]}
