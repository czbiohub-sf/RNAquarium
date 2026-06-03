#!/bin/bash

#SBATCH --job-name=stepv7b_RNAquarium_metatranscriptome_viralclustering_AFTERBBdukfilterandtrim
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
## taxonomy_hits_viruses_nonphage_mostrecent.tsv


#module purge
#module load anaconda/2023.03
#conda activate seqtk
#time_stamp=$(date +%Y_%m_%d)

#mv RNAquarium_outputs/list_of_targetids_forclustering.txt list_of_targetids_forclustering.txt
#mv RNAquarium_outputs/allchunks_blastnanddiamond_hits_viruses_mostrecent.fasta allchunks_blastnanddiamond_hits_viruses_mostrecent.fasta


cd RNAquarium_outputs/virus_outputs
date
efetch -input "list_of_targetids_forclustering.txt" -db nuccore -format fasta >> targetNT_sequences0.fasta # see below we want to delete this before next run
sleep 1
date
## we want to replace all spaces, commas, etc with underscores
sed 's/[^a-zA-Z0-9>]/_/g' targetNT_sequences0.fasta > targetNT_sequences1.fasta
sed 's/__/_/g' targetNT_sequences1.fasta > targetNT_sequences.fasta
sleep 1
date ## below not nonphage, but all!
cat taxonomy_hits_viruses_mostrecent.fasta targetNT_sequences.fasta > taxonomy_hits_viruses_mostrecent_withtargets.fasta
rm targetNT_sequences0.fasta
date
sleep 1


## IF USING BBDUK FILTER & STEP 1B, THEN ALSO RUN STEP 7B HERE TO FILTER & TRIM ALL VIRAL SEQUENCES BEFORE CLUSTERING
## saving initial fasta as _prefilter.fasta
cp taxonomy_hits_viruses_mostrecent_withtargets.fasta taxonomy_hits_viruses_mostrecent_withtargets_prefilter.fasta

## first BBduk filter of adapters/primers etc.

module load anaconda
conda activate bbmap

#cd RNAquarium_outputs/virus_outputs # already here!
# Copy adapter/primer list to this directory
cp /path/to/RNAquarium/src/metatranscriptome/bin/fastp_adapters_with9added.fasta .

## running filter only on viruses, not NTtargets
bbduk.sh -Xmx24g threads=auto in=taxonomy_hits_viruses_mostrecent.fasta out=taxonomy_hits_viruses_mostrecent_masked.fa ref=adapters,artifacts,phix,lambda,pjet,kapa,fastp_adapters_with9added.fasta k=21 mink=10 hdist=1 maskmiddle=f kmask=N fastawrap=0 stats=taxonomy_hits_viruses_mostrecent_masking_stats.txt

# Step 2: Remove all Ns
sed '/^>/!s/N//g' taxonomy_hits_viruses_mostrecent_masked.fa > taxonomy_hits_viruses_mostrecent_masked_trimmed.fa 

## now combine filtered virus fasta with NTtargets
cat taxonomy_hits_viruses_mostrecent_masked_trimmed.fa targetNT_sequences.fasta > taxonomy_hits_viruses_mostrecent_withtargets.fasta

sleep 1
### now mmseqs cluster

mmseqs easy-cluster taxonomy_hits_viruses_mostrecent_withtargets.fasta taxonomy_hits_viruses_easycluster1 tmp --min-seq-id 0.5 -c 0.1 --cov-mode 1



#mv list_of_targetids_forclustering.txt RNAquarium_outputs/list_of_targetids_forclustering.txt
#mv targetNT_sequences.fasta RNAquarium_outputs/targetNT_sequences.fasta
#mv allchunks_blastnanddiamond_hits_viruses_mostrecent.fasta RNAquarium_outputs/allchunks_blastnanddiamond_hits_viruses_mostrecent.fasta

#mv RNAquarium_outputs/list_of_targetids_forclustering.txt list_of_targetids_forclustering.txt
