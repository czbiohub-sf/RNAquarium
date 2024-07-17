#!/bin/bash
#SBATCH --job-name=zf-nonhost-nf
#SBATCH --time=14-12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -e slurm.out/slurm-%A_%a.err
#SBATCH -o slurm.out/slurm-%A_%a.out

module load nextflow
module load mamba

set -o verbose
NXF_OPTS="-Xms500M -Xmx16G" PATH=$PATH:$PWD/bin nextflow run main.nf \
	-params-file params.example.yaml \
	-profile bruno,mamba