#!/bin/bash

#SBATCH --job-name=step3_RNAquarium_metatranscriptome_addsequencetotsv
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1300G
#SBATCH --partition cpu         # Partition to submit
#SBATCH --gpus 0            # Reserve 0 GPUs for usage
#SBATCH --cpus-per-task=32
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-user your.email@example.com   # this is the email you wish to be notified at
#SBATCH --mail-type END,FAIL   # ALL will alert you of job beginning, completion, failure etc

## just the bash commands I would use

module purge
module load R/4.4
#R
#source("nr_diamond_processing.r")
Rscript 3_addsequence.r