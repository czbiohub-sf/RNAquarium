#!/bin/bash

#SBATCH --job-name=step2_RNAquarium_metatranscriptome_createfastas
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=300G
#SBATCH --partition cpu         # Partition to submit
#SBATCH --gpus 0            # Reserve 0 GPUs for usage
#SBATCH --cpus-per-task=16
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-user eric.waltari@czbiohub.org   # this is the email you wish to be notified at
#SBATCH --mail-type END,FAIL   # ALL will alert you of job beginning, completion, failure etc

## just the bash commands I would use

module purge
module load anaconda/2023.03
conda activate seqtk

mv RNAquarium_outputs/*list.txt .
## nonzfhum.fasta chunk_*_nonzfhum.fasta
for fasta in chunk_*_nonzfhum.fasta; do
  date
  seqtk subseq "$fasta" taxonomy_hits_nonhost_list.txt >> taxonomy_hits_nonhost_list.fasta
done

mv *list.fasta RNAquarium_outputs/
mv *list.txt RNAquarium_outputs/