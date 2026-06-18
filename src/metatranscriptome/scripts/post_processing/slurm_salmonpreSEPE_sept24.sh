#!/bin/bash

#SBATCH --job-name=salmon_setup_sept24
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=235G
#SBATCH --partition cpu         # Partition to submit
#SBATCH --gpus 0            # Reserve 0 GPUs for usage
#SBATCH --cpus-per-task=16
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-user your.email@example.com   # this is the email you wish to be notified at
#SBATCH --mail-type END,FAIL   # ALL will alert you of job beginning, completion, failure etc

module purge
module load anaconda/2023.03
conda activate salmon
cd /path/to/salmon_steps/nonhost_viruscounts_salmon

salmon index -t taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-17_sprivireplacediwth2genomes_masked_trimmed.fa -i salmon_index_sept24 -k 31 --keepDuplicates
