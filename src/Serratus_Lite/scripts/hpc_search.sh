#!/bin/bash

#SBATCH --job-name=Serratus
#SBATCH --time=14-00:00:00
#SBATCH --array=1-54189%100
#SBATCH --partition preempted
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -e slurm.out/slurm-%j.err
#SBATCH -o slurm.out/slurm-%j.out

basedir="${1}"
accessions_list="${2}"
seqbasedir="${3}"

virus_files="${basedir}/query_virus_data/*.fa"

# declare arrays
readarray -t accessions < <(cat "${accessions_list}")

declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))

# #working dir
# working_dir=${basedir}/Serratus_Lite_pipeline/zpoxv # <---- adjust this when switching query!!!!
# cd $working_dir

# #setup query
# query=zpoxv.fasta # <---- adjust this when switching query!!!!

# run bowtie2
# -f <input is fasta>  -U <input is unpaired>
for VIRALSEQ in ${virus_files}
do
    query="$(basename "${VIRALSEQ}")"
    virus_name="$(basename "${VIRALSEQ%.*}")"
    working_dir="${basedir}/results/${virus_name}"
    FQ1="${seqbasedir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq.gz"
    FQ2="${seqbasedir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.gsnapFiltered.fastq.gz"

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
    # SUMZER_COMMENT=$(echo sra="na",genome="${query}",version=,date=$(date +%y%m%d-%R))
    SUMZER_TSV="${basedir}/query_virus_data/${query}.sumzer.tsv"
    summarizer="python3 serratus_summarizer.py /dev/stdin ${SUMZER_TSV} $working_dir/summarized/${accessions[$idx]}.summary /dev/null"

    if [ -e $working_dir/bam/${accessions[$idx]}.bam ]
    then
        samtools view $working_dir/bam/${accessions[$idx]}.bam | $summarizer
        rm $working_dir/bam/${accessions[$idx]}.bam
    fi
done




