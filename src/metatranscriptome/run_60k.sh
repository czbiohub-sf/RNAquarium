#!/usr/bin/env bash
#SBATCH --job-name=nf-RNAquarium
#SBATCH --chdir=/hpc/mydata/gibraan.rahman/projects/RNAquarium/merge_assemble
#SBATCH --output=./slurm/nf-RNAquarium.out
#SBATCH --time=21-00:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=cpu
#SBATCH --mail-user=gibraan.rahman@czbiohub.org
#SBATCH --mail-type=BEGIN,END,FAIL

UNMERGED_ACC="/hpc/projects/balla_group/sra_experiments/zebrafish_RNAseq/rnaquarium_output/nonhost_reads/"
BIOPROJ_MAP="tmp/filt_mapping.json"
WORK_DIR="/hpc/scratch/group.swe/rnaquarium/work"
PUB_DIR="./output_2024_07_25/"
NT_DIR="/hpc/scratch/group.swe/db/nt_clustered"
NT_DB_NAME="nt_compressed_shuffled.fa"
NTFULL_DIR="/hpc/scratch/group.data.science/eric_temp/databases/2025/core_nt"
NTFULL_DB_NAME="core_nt"
NR_DIR="/hpc/scratch/group.data.science/eric_temp/databases/2025/nr_clustered"
TAXONOMY_DB="/hpc/scratch/group.data.science/eric_temp/mmseqs_out/taxonomizr/feb2025taxonomy/nameNode.sqlite"

mkdir -p $PUB_DIR

REPORT="${PUB_DIR}/report.60k.html"
TIMELINE="${PUB_DIR}/timeline.60k.html"

nextflow run \
    -profile slurm,singularity \
    -with-report $REPORT \
    -with-timeline $TIMELINE \
    -w $WORK_DIR \
    -ansi-log false \
    -resume 76dec62b-0cac-4447-9818-912069d370a9 \
    --bioproj_map $BIOPROJ_MAP \
    --unmerged_accessions $UNMERGED_ACC \
    --publish_dir $PUB_DIR \
    --nt_dir $NT_DIR \
    --nt_db_name $NT_DB_NAME \
    --nr_dir $NR_DIR \
    --taxonomy_db $TAXONOMY_DB \
    main.nf
