#!/usr/bin/env bash
#SBATCH --job-name=nf-RNAquarium
#SBATCH --chdir=/path/to/RNAquarium/src/metatranscriptome
#SBATCH --output=./nf-RNAquarium.out
#SBATCH --time=21-00:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=cpu
#SBATCH --mail-user=your.email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL

UNMERGED_ACC="/path/to/RNAquarium/src/metatranscriptome/all_unmapped_links/"
BIOPROJ_MAP="/path/to/SRA_metadata/bioproject_mapping.json"
WORK_DIR="/path/to/scratch/nf_metatmp/"
PUB_DIR="/path/to/output/metatranscriptome/"
NT_DIR="/path/to/databases/nt_clustered"
NT_DB_NAME="nt_compressed_shuffled.fa"
NTFULL_DIR="/path/to/databases/2025/core_nt"
NTFULL_DB_NAME="core_nt"
NR_DIR="/path/to/databases/2025/nr_clustered"
TAXONOMY_DB="/path/to/databases/taxonomizr/feb2025taxonomy/nameNode.sqlite"

mkdir -p $PUB_DIR

REPORT="${PUB_DIR}/report.75k.html"
TIMELINE="${PUB_DIR}/timeline.75k.html"

unset NXF_VER

module load nextflow/24.10.5
## if module load is unavailable, you can call the nextflow binary by full path; note the unset NXF_VER above
nextflow run \
    -profile slurm,singularity \
    -with-report $REPORT \
    -with-timeline $TIMELINE \
    -w $WORK_DIR \
    -ansi-log false \
    --bioproj_map $BIOPROJ_MAP \
    --unmerged_accessions $UNMERGED_ACC \
    --publish_dir $PUB_DIR \
    --nt_dir $NT_DIR \
    --nt_db_name $NT_DB_NAME \
    --ntfull_dir $NTFULL_DIR \
    --ntfull_db_name $NTFULL_DB_NAME \
    --nr_dir $NR_DIR \
    --taxonomy_db $TAXONOMY_DB \
    -resume \
    main.nf
