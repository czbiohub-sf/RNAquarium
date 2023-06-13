#purge output folders
rm -rf slurm.out
mkdir -p slurm.out
rm -rf summarized
mkdir -p summarized
rm -rf bam
mkdir -p bam

#submit slurm job array
sbatch search_nonZFv2.sh