#!/bin/bash

#SBATCH --job-name=stepv5A_RNAquarium_metatranscriptome_viruscuration_justrealmfix
#SBATCH --time=21:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=3900G
#SBATCH --partition cpu         # Partition to submit
#SBATCH --gpus 0            # Reserve 0 GPUs for usage
#SBATCH --cpus-per-task=8
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-user eric.waltari@czbiohub.org   # this is the email you wish to be notified at
#SBATCH --mail-type END,FAIL   # ALL will alert you of job beginning, completion, failure etc

## just the bash commands I would use

module purge
module load R/4.4
#R
#source("nr_diamond_processing.r")
Rscript virus_curation_plots_addsequenceA.r

## this changes names of files for running virus 1b if needed
sleep 2
cd RNAquarium_outputs
mv taxonomy_hits_viruses0_fullcols_mostrecent.tsv taxonomy_hits_viruses0_fullcols_mostrecent_norealmcolumn.tsv
cd virus_outputs
mv taxonomy_hits_viruses0_fullcols_mostrecent2.tsv ../taxonomy_hits_viruses0_fullcols_mostrecent.tsv
