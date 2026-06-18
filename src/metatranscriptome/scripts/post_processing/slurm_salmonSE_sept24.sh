#!/bin/bash

#SBATCH --job-name=salmon_SEportion_sept24
#SBATCH --time=2-21:00:00      # NOTE THAT 22 HOURS IS A CUT-OFF, USING LESS HERE
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
log_file="salmonlogSE_${current_date}.txt"

# Loop through all subfolders directly inside "Paired"
## now using actual location of unmapped_reads since we are no longer working in that folder
for folder in /path/to/host_mapping/unmapped_reads/Single/*/; do
  # Find the fastq file in the current folder
  fastq=$(find "$folder" -type f -name "*.fastq.gz")

  # Extract the subfolder name
  subfolder_name=$(basename "$folder")

  # Check if a FASTQ file was found
  if [[ -n "$fastq" ]]; then
    # Run Salmon quant with single-end input
    $RUN salmon quant \
      -i /path/to/salmon_steps/nonhost_viruscounts_salmon/salmon_index_sept24 \
      -l U \
      -r ${fastq} \
      -p 32 \
      -o salmon_transcripts_sept24/${subfolder_name}_salmon \
      --validateMappings

    echo "✅ Processed ${fastq} → Output: salmon_transcripts_sept24/${subfolder_name}_salmon" | tee -a $log_file
  else
    echo "⚠️  No FASTQ file found in ${folder}" | tee -a $log_file
  fi
done
