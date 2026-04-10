#!/bin/bash

#SBATCH --job-name=step2b_RNAquarium_metatranscriptome_bbkukmaskingfilterandflag
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --partition cpu         # Partition to submit
#SBATCH --gpus 0            # Reserve 0 GPUs for usage
#SBATCH --cpus-per-task=16
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-user eric.waltari@czbiohub.org   # this is the email you wish to be notified at
#SBATCH --mail-type END,FAIL   # ALL will alert you of job beginning, completion, failure etc

## just the bash commands I would use

module purge
module load anaconda
conda activate bbmap

cd RNAquarium_outputs

# Copy adapter/primer list to this directory
cp /hpc/projects/balla_group/sra_experiments/RNAquarium_75k/RNAquarium/src/metatranscriptome/fastp_adapters_with9added.fasta .

# Check if masked.fasta exists, only continue if not & also save previous version of unmasked fasta!
if [ ! -f taxonomy_hits_nonhost_list_masked.fasta ]; then
    cp taxonomy_hits_nonhost_list.fasta taxonomy_hits_nonhost_list_unmasked.fasta
else
    echo "BBDuk was already run on this version of taxonomy_hits_nonhost_list.fasta"
    exit 1
fi


bbduk.sh -Xmx16g threads=auto in=taxonomy_hits_nonhost_list.fasta out=taxonomy_hits_nonhost_list_masked.fasta ref=adapters,artifacts,phix,lambda,pjet,kapa,fastp_adapters_with9added.fasta k=21 mink=10 hdist=1 maskmiddle=f kmask=N fastawrap=0 statscolumns=5 stats=taxonomy_hits_nonhost_masking_stats.txt
#This again uses both BBDuk set of adapters and artifacts but also the expanded list from Fastp set of adapters + the additional I have found (e.g. Clontech primers)

sleep 9
#Next is the command to flag all sequences that have less than 150 non-N bases.
awk '/^>/ {name=$0; next} {gsub(/N/,""); if(length($0)<150) print name}' taxonomy_hits_nonhost_list_masked.fasta | sed 's/^>//' > taxonomy_hits_nonhost_sequences_to_exclude.txt
sleep 2
