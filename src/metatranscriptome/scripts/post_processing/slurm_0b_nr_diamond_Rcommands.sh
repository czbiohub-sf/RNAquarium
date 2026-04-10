#!/bin/bash

#SBATCH --job-name=step0b_RNAquarium_metatranscriptome_rawdiamondNRprocessing_main
#SBATCH --time=3-22:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1300G
#SBATCH --partition cpu         # Partition to submit
#SBATCH --gpus 0            # Reserve 0 GPUs for usage
#SBATCH --cpus-per-task=32
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-user eric.waltari@czbiohub.org   # this is the email you wish to be notified at
#SBATCH --mail-type END,FAIL   # ALL will alert you of job beginning, completion, failure etc

## just the bash commands I would use

cd nr_diamond
module purge
module load R/4.4
#R
#source("nr_diamond_processing.r")
Rscript nr_diamond_processing.r


# Concatenate all fishorprimates files into one combined file
cat *_fishorprimates.txt > chunks_diamond_hits_fishorprimatescombined.txt

# Verify the combined file was created successfully before deleting intermediates
if [ -s chunks_diamond_hits_fishorprimatescombined.txt ]; then
    echo "Combined file created successfully with $(wc -l < chunks_diamond_hits_fishorprimatescombined.txt) lines"
    # Delete the intermediate files
    rm *_fishorprimates.txt
    echo "Intermediate fishorprimates files deleted"
else
    echo "ERROR: Combined file is empty or was not created. Intermediate files NOT deleted."
    exit 1
fi