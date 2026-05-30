#!/bin/bash

#SBATCH --job-name=diamond_taxid_update_usingblastdbcmd_full216diamondruns
#SBATCH --time=22:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --partition cpu         # Partition to submit
#SBATCH --gpus 0            # Reserve 0 GPUs for usage
#SBATCH --cpus-per-task=4
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-user your.email@example.com   # this is the email you wish to be notified at
#SBATCH --mail-type END,FAIL   # ALL will alert you of job beginning, completion, failure etc

## just the bash commands I would use

## move to nr_diamond folder?? no just run from nr_diamond
#cd nr_diamond

echo "Starting diamond_taxid_update_blastdbcmd.py script at $(date)"
echo "This script updates Diamond search results with missing taxids using your local BlastDBCmd"


module load anaconda
python -u diamond_taxid_update_blastdbcmd.py

echo "Starting diamond_taxid_moveold_andrename.py script at $(date)"
echo "This script moves and renames the updated Diamond result files to efetch_update/ with original .diamond.txt.gz names"

module load anaconda
python -u diamond_taxid_moveold_andrename.py

## now move to nr_diamond folder and see if there are .diamond.txt.gz files in efetch_update folder
## if there are, first move all .diamond.txt.gz files in this folder to new folder called month_date/
## then move all .diamond.txt.gz files from the efetch_update folder to this folder

# Move into the nr_diamond directory
#cd nr_diamond || { echo "ERROR: Failed to enter nr_diamond directory"; exit 1; }

# Check if any .diamond.txt.gz files exist in efetch_update/


# Check for .diamond.txt.gz files in efetch_update/
if ls efetch_update/*.diamond.txt.gz 1> /dev/null 2>&1; then
  # Step 1: Generate folder name based on month and year, e.g., Oct_2025
  BASE_NAME=$(date +%b_%Y)
  FOLDER_NAME="$BASE_NAME"
  COUNT=2

  # Step 2: Ensure the folder name is unique by appending _2, _3, etc.
  while [ -d "$FOLDER_NAME" ]; do
    FOLDER_NAME="${BASE_NAME}_$COUNT"
    COUNT=$((COUNT + 1))
  done

  # Step 3: Create the unique folder
  mkdir "$FOLDER_NAME"
  echo "Created new folder: $FOLDER_NAME"

  # Step 4: Move any *.diamond.txt.gz files already in nr_diamond/ into the new folder
  mv -v ./*.diamond.txt.gz "$FOLDER_NAME/" 2>/dev/null
  echo "Moved existing .diamond.txt.gz files into $FOLDER_NAME/"

  # Step 5: Move *.diamond.txt.gz files from efetch_update/ up into nr_diamond/
  mv -v efetch_update/*.diamond.txt.gz ./ 2>/dev/null
  echo "Pulled new .diamond.txt.gz files from efetch_update/ into nr_diamond/"

else
  echo "No .diamond.txt.gz files found in efetch_update/"
fi

