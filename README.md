# RNAquarium merge & assembly steps

## Example usage

Test Data:

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

Using conda instead of singularity

```
nextflow run -profile test,slurm,conda main.nf
```
