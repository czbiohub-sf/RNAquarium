#!/bin/bash

#SBATCH --job-name=salmon_PEportion_sept24
#SBATCH --time=2-22:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1260G
#SBATCH --partition cpu         # Partition to submit
#SBATCH --gpus 0            # Reserve 0 GPUs for usage
#SBATCH --cpus-per-task=32
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-user your.email@example.com   # this is the email you wish to be notified at
#SBATCH --mail-type END,FAIL   # ALL will alert you of job beginning, completion, failure etc

## just the bash commands I would use srun --pty --cpus-per-task 28 --mem=1260G --time=1-10 bash -l

#module purge
module purge
module load anaconda/2023.03
conda activate salmon

cd /path/to/salmon_steps/

#cd /path/to/host_mapping/unmapped_reads/

## new location of indexes
# /path/to/salmon_steps/nonhost_viruscounts_salmon

# Get the current date and format it (e.g., YYYY-MM-DD)
current_date=$(date +%Y-%m-%d)

# Define the log file name with the date included
log_file="salmonlogPE_${current_date}.txt"

# Loop through all immediate subfolders inside "Paired"
for folder in /path/to/host_mapping/unmapped_reads/Paired/*/; do
  # Find the mate1 and mate2 fastq files in the current folder
  mate1=$(find "$folder" -type f -name "Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq.gz")
  mate2=$(find "$folder" -type f -name "Unmapped.out.mate2.filteredbyBT.dedup.gsnapFiltered.fastq.gz")

  # Extract the subfolder name
  subfolder_name=$(basename "$folder")

  # Check if both mate1 and mate2 files exist
  if [[ -n "$mate1" && -n "$mate2" ]]; then
    # Run the salmon quant command with both mates
    $RUN salmon quant -i /path/to/salmon_steps/nonhost_viruscounts_salmon/salmon_index_sept24 -l IU \
        -1 ${mate1} \
        -2 ${mate2} \
        -p 32 \
        -o salmon_transcripts_sept24/${subfolder_name}_salmon \
        --validateMappings

    echo "✅ Processed ${mate1} and ${mate2} → Output: salmon_transcripts_sept24/${subfolder_name}_salmon" | tee -a $log_file
  else
    echo "⚠️  Mate1 or Mate2 fastq files missing in ${folder}" | tee -a $log_file
  fi
done

