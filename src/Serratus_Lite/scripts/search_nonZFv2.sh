#!/bin/bash

#SBATCH --job-name=Zpoxv_Serratus
#SBATCH --time=14-00:00:00
#SBATCH --array=1-44667%1000
#SBATCH --partition preempted
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -e slurm.out/slurm-%j.err
#SBATCH -o slurm.out/slurm-%j.out

# declare arrays
readarray -t accessions < <(cat /hpc/projects/balla_group/sra_experiments/all_zebrafish_RNAseq/SRA_accession_list.txt)

declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))

module load anaconda
conda activate bowtie2

#setting directories
basedir=/hpc/projects/balla_group/sra_experiments/search_nonZFv2
seqbasedir=/hpc/projects/balla_group/sra_experiments/all_zebrafish_RNAseq/unmapped_dev/Gsnap_out
#working dir
working_dir=${basedir}/Serratus_Lite_pipeline/zpoxv # <---- adjust this when switching query!!!!
cd $working_dir

#setup query
query=zpoxv.fasta # <---- adjust this when switching query!!!!

# run bowtie2
# -f <input is fasta>  -U <input is unpaired>

FQ1=${seqbasedir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq.gz
FQ2=${seqbasedir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.gsnapFiltered.fastq.gz
VIRALSEQ="${basedir}/query_virus_data/${query}"

if [ -e $FQ1 ] && [ -e $FQ2 ]
then
    #PE
    bowtie2 --quiet --very-sensitive-local \
    --rg-id na --rg LB:na --rg SM:na \
    --rg PL:na --rg PU:na \
    -x ${VIRALSEQ} -1 $FQ1 -2 $FQ2 -q | \
    samtools view -b -F 4 - > $working_dir/bam/${accessions[$idx]}.bam
else
    if [ -e $FQ1 ]
    then
        #SE
        bowtie2 --quiet --very-sensitive-local \
        --rg-id na --rg LB:na --rg SM:na \
        --rg PL:na --rg PU:na \
        -x ${VIRALSEQ} -U $FQ1 -q | \
        samtools view -b -F 4 - > $working_dir/bam/${accessions[$idx]}.bam
    fi
fi

#Summarizer
SUMZER_COMMENT=$(echo sra="na",genome="${query}",version=,date=$(date +%y%m%d-%R))
SUMZER_TSV="${basedir}/query_virus_data/${query}.sumzer.tsv"
summarizer="python3 serratus_summarizer.py /dev/stdin ${SUMZER_TSV} $working_dir/summarized/${accessions[$idx]}.summary /dev/null"

if [ -e $working_dir/bam/${accessions[$idx]}.bam ]
then
    samtools view $working_dir/bam/${accessions[$idx]}.bam | $summarizer
    rm $working_dir/bam/${accessions[$idx]}.bam
fi





