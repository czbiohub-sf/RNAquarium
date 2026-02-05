#!/usr/bin/env bash
#SBATCH --job-name=nf-RNAquarium
#SBATCH --chdir=/hpc/projects/balla_group/sra_experiments/RNAquarium_75k/RNAquarium/src/metatranscriptome
#SBATCH --output=./nf-RNAquarium.out
#SBATCH --time=21-00:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=cpu
#SBATCH --mail-user=eric.waltari@czbiohub.org
#SBATCH --mail-type=BEGIN,END,FAIL

UNMERGED_ACC="/hpc/projects/balla_group/sra_experiments/RNAquarium_75k/RNAquarium/src/metatranscriptome/all_unmapped_links/"
BIOPROJ_MAP="/hpc/projects/balla_group/sra_experiments/SRA_metadata/bioproject_mapping.json"
WORK_DIR="/hpc/scratch/group.data.science/eric_temp/nf_metatmp_24jul25/"
PUB_DIR="/hpc/projects/balla_group/sra_experiments/versioned_zf_output/75k_unstable/metatranscriptome/"
NT_DIR="/hpc/scratch/group.swe/db/nt_clustered"
NT_DB_NAME="nt_compressed_shuffled.fa"
NTFULL_DIR="/hpc/scratch/group.data.science/eric_temp/databases/2025/core_nt"
NTFULL_DB_NAME="core_nt"
NR_DIR="/hpc/scratch/group.data.science/eric_temp/databases/2025/nr_clustered"
TAXONOMY_DB="/hpc/scratch/group.data.science/eric_temp/mmseqs_out/taxonomizr/feb2025taxonomy/nameNode.sqlite"

mkdir -p $PUB_DIR

REPORT="${PUB_DIR}/report.75k.html"
TIMELINE="${PUB_DIR}/timeline.75k.html"

unset NXF_VER

module load nextflow/24.10.5
## to run this version /hpc/apps/x86_64/nextflow/24.10.5/bin/nextflow /path/to/module/nextflow/24.10.5/bin/nextflow run but getting errors due to my .bash_profile note new unset command first
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
