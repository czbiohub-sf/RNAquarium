#!/bin/bash

#purge output folders
virus_results="/hpc/scratch/group.theory/jparas/Zebrafish-RNA-Quarium/src/Serratus_Lite/results/${1}"
mkdir $virus_results
cd $virus_results #TODO rel path
rm -rf slurm.out
mkdir -p slurm.out
rm -rf summarized
mkdir -p summarized
rm -rf bam
mkdir -p bam

#submit slurm job array
sbatch search_nonZFv2.sh $1