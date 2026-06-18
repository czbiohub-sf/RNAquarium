#!/bin/bash

#SBATCH --job-name=salmon_outputs_collating_SEPE_sept24
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=260G
#SBATCH --partition cpu         # Partition to submit
#SBATCH --gpus 0            # Reserve 0 GPUs for usage
#SBATCH --cpus-per-task=28
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-user your.email@example.com   # this is the email you wish to be notified at
#SBATCH --mail-type END,FAIL   # ALL will alert you of job beginning, completion, failure etc


module purge
module load anaconda/2023.03
conda activate salmon
cd /path/to/salmon_steps


#Collating Salmon outputs step
#cd /path/to/rnaquarium_output/nonhost_reads_split/single

## note this will reset all.quant0.tsv each time, and then all.quantt_mostrecent.tsv each time
## this means file needs to be renamed at end with timestamp
first_file=$(find salmon_transcripts_sept24/SRR102*_salmon/quant.sf | head -n 1)
cut -f1 ${first_file} | tr '\n' '\t' > all.quant0.tsv
echo "" >> all.quant0.tsv
sed 's/Name/run_name/' all.quant0.tsv > all.quantt_sept24.tsv

for file in salmon_transcripts_sept24/*_salmon/quant.sf; do
  # Extract the folder name
  folder_name=$(basename $(dirname $file))
  
  # Print the folder name as the first column in the new row
  printf "%s\t" ${folder_name} >> all.quantt_sept24.tsv
  
  # Extract NumReads column, convert it to a row (tab-separated), and append to all.quantt
  tail -n+2 ${file} | cut -f5 | tr '\n' '\t' >> all.quantt_sept24.tsv
  echo "" >> all.quantt_sept24.tsv  # Add a newline at the end of each row
done

gzip all.quantt_sept24.tsv

# Final Salmon count matrix (~75k SRA runs x ~182k virus contigs):
#   <OUTPUT_ROOT>/host_mapping/unmapped_reads/all.quantt_sept24.tsv.gz

echo "Transposed summary saved to all.quantt_sept24.tsv"



