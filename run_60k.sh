UNMERGED_ACC="/hpc/projects/balla_group/sra_experiments/RNAquarium_prototyping/rnaquarium_unmapped_sample"
BIOPROJ_MAP="tmp/full_mapping.json"
WORK_DIR="/hpc/scratch/group.swe/rnaquarium/work"

nextflow run \
    -profile slurm,singularity \
    -with-report report.60k.html \
    -with-timeline timeline.60k.html \
    -w $WORK_DIR \
    -resume \
    --bioproj_map $BIOPROJ_MAP \
    --unmerged_accessions  $UNMERGED_ACC \
    main.nf
