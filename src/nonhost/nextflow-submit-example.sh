#!/bin/bash
#SBATCH --job-name=zf-nonhost-nf
#SBATCH --time=14-12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -e slurm.out/slurm-%A_%a.err
#SBATCH -o slurm.out/slurm-%A_%a.out

module load nextflow
module load mamba

ACCESSIONS=data/SRA_accessions_list.test.txt

set -o verbose
NXF_OPTS="-Xms500M -Xmx16G" PATH=$PATH:$PWD/bin nextflow run main.nf \
	--accession-list $ACCESSIONS \
        --ref-genome /hpc/projects/balla_group/sra_experiments/tools/STAR/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa \
        --ref-genome-gtf /hpc/projects/balla_group/sra_experiments/tools/STAR/Danio_rerio.GRCz11.108.gtf \
        --ercc-fa /hpc/projects/balla_group/sra_experiments/tools/STAR/ERCC92.fa \
        --ercc-gtf /hpc/projects/balla_group/sra_experiments/tools/STAR/ERCC92.gtf \
	--genome-size 1396431182 --skip-host-counts false \
	--hisat-use-transcript true \
	--star-threads-small 4 --star-threads-large 16 --max-cpus 128 --max-memory "256GB" \
	--star-use-shared-mem false --parallel-downloads 50 \
	--tmp /tmp/ \
	--cleanup-intermediate true --nxf-unstage-hack true \
	-profile bruno,mamba


#	--fastq-path $PWD/fastq \
# --backup-tmp /hpc/scratch/<your_temp> --backupScratchHack true
