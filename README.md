# RNAquarium merge & assembly steps

Example usage with test data:

```
nextflow run -profile test,bruno main.nf
```

Example usage with subset of test data (only two accessions):

```
nextflow run -profile test,bruno \
    --accession_list test/accession_subset.txt \
    --bioproj_map false \
    main.nf
```
