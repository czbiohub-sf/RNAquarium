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
readarray -t accessions < <(cat "${ACCESSIONS_LIST}") #54189 64G 4cpu

#output directories
fdir=${WORKING_DIR}/fastq
pdir=${WORKING_DIR}/prefetched
sdir=${WORKING_DIR}/STAR_out

echo "accession: ${accessions[$idx]}"

############################
# proceed with current job #
############################
echo "processing..." >> ${WORKING_DIR}/logs/STAR.process.log
#sleep 0-120 seconds to space out 100 concurrent jobs from doing prefetch at the same time (prevent timing outs)
sleep $((idx % 120))
cd $WORKING_DIR

# make directories for non STAR stdout and stderr.  STAR out directories for this accession will be purged later in the STAR section
if [ -e "${sdir}/other_stdout_stderr/${accessions[$idx]}" ]
then
    rm -rf ${sdir}/other_stdout_stderr/${accessions[$idx]}
fi
mkdir ${sdir}/other_stdout_stderr/${accessions[$idx]}

###############
# fasterqdump #
###############
#create folder for fastq
if [ -e ${fdir}/${accessions[$idx]} ]
then
    rm -rf ${fdir}/${accessions[$idx]}
fi
mkdir ${fdir}/${accessions[$idx]}

#prefetch
${SRA_BIN}/prefetch --max-size 1t --force ALL --output-directory ${pdir} ${accessions[$idx]} 1> ${sdir}/other_stdout_stderr/${accessions[$idx]}/prefetch.stdout.txt 2> ${sdir}/other_stdout_stderr/${accessions[$idx]}/prefetch.stderr.txt

#fasterq-dump
cd ${pdir}
${SRA_BIN}/fasterq-dump --split-3 --mem 24G --outdir ${fdir}/${accessions[$idx]} ${accessions[$idx]} 1> ${sdir}/other_stdout_stderr/${accessions[$idx]}/fqdump.stdout.txt 2> ${sdir}/other_stdout_stderr/${accessions[$idx]}/fqdump.stderr.txt  && {
    echo fasterq-dump: no error
} || {
    echo fasterq-dump encounterd error, reverting back to using fastqdump
    ${SRA_BIN}/fastq-dump --split-3 --disable-multithreading --outdir ${fdir}/${accessions[$idx]} ${accessions[$idx]} 1> ${sdir}/other_stdout_stderr/${accessions[$idx]}/fqdump.stdout.txt 2> ${sdir}/other_stdout_stderr/${accessions[$idx]}/fqdump.stderr.txt
}

#remove prefretched data
cd ..
rm -rf ${pdir}/${accessions[$idx]}
