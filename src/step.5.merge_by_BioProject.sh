#!/bin/bash

#SBATCH --job-name=merge_fq
#SBATCH --time=14-00:00:00
#SBATCH --array=1-1498%500
#SBATCH --nodes=1
#SBATCH --partition preempted
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH -e slurm.out.merge_fq_by_BioProject/slurm-%A_%a.err
#SBATCH -o slurm.out.merge_fq_by_BioProject/slurm-%A_%a.out

declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))

module load anaconda
conda activate zf_pipeline

#setting directories
#working_dir="/hpc/projects/balla_group/sra_experiments/all_zebrafish_RNAseq/unmapped_dev"
working_dir="/hpc/scratch/group.theory/jparas/zf_pipeline"

gdir=${working_dir}/Gsnap_out
outdir=${working_dir}/out_fq_by_BioProjects
key=BioProject

# declare arrays
readarray -t IDs < <(cat "${working_dir}/SRA_accession_list.1.27.23.txt_BioProjects.txt")

#skip if results exists
# if [ -e ${outdir}/${IDs[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq.gz ]
# then
#     echo "skipping ${IDs[$idx]} because merged fastq.gz file exists"
#     exit 0
# fi

#clear previous results
if [ -d ${outdir}/${IDs[$idx]} ]
then
    rm -rf ${outdir}/${IDs[$idx]}
fi

#main command
python merge_fq_by_ID.py --ID ${IDs[$idx]} --key ${key}  --indir ${gdir} --outdir ${outdir}
