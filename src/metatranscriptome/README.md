# RNAquarium merge & assembly steps

## Setup

Make sure you have built the Apptainer image required for running the taxonomy scripts.

```
apptainer build \
    src/metatranscriptome/containers/taxonomy.sif \
    src/metatranscriptome/containers/taxonomy.def
```

## Example usage

Test data:

```
nextflow run -profile test,slurm,singularity main.nf
```

Subset of test data (only two accessions):

```
nextflow run -profile test,slurm,singularity \
    --accession_list test/accession_subset.txt \
    --bioproj_map false \
    main.nf
```

Using conda instead of singularity:

```
nextflow run -profile test,slurm,conda main.nf
```
