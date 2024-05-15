#!/usr/bin/env bash
#SBATCH --job-name=nf-RNAquarium
#SBATCH --chdir=/hpc/mydata/gibraan.rahman/projects/RNAquarium/merge_assemble
#SBATCH --output=./slurm/nf-RNAquarium.out
#SBATCH --time=3-00:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=cpu
#SBATCH --mail-user=gibraan.rahman@czbiohub.org
#SBATCH --mail-type=BEGIN,END,FAIL

UNMERGED_ACC="/hpc/projects/balla_group/sra_experiments/zebrafish_RNAseq/rnaquarium_output/nonhost_reads/"
BIOPROJ_MAP="tmp/full_mapping.json"
WORK_DIR="/hpc/scratch/group.swe/rnaquarium/work"
PUB_DIR="./output/"

mkdir -p $PUB_DIR

REPORT="${PUB_DIR}/report.60k.html"
TIMELINE="${PUB_DIR}/timeline.60k.html"

nextflow run \
    -profile slurm,singularity \
    -with-report $REPORT \
    -with-timeline $TIMELINE \
    -w $WORK_DIR \
    -ansi-log false \
    -resume 97b8256e-2675-40da-9f81-23c31ddbc327 \
    -process.cache=lenient \
    --bioproj_map $BIOPROJ_MAP \
    --unmerged_accessions $UNMERGED_ACC \
    --publish_dir $PUB_DIR \
    main.nf
