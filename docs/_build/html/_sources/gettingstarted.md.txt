# Getting Started
## Dependencies
The RNAquarium pipeline requires several tools to be available in `$PATH`.

They can either be installed separately, or, more typically, provided by containers.
`-profile docker` or `-profile singularity`

If using containers, there are a small number of programs that are not available as containers.
These can be installed with
```bash
nextflow run install-deps.nf
```

## Running the pipeline
the main preprocessing pipeline is the default entry point of `main.nf`.

The minimum required inputs are 
 - a host genome
 - host genome annotations
 - an estimate of the host genome size in bytes
 - ERCC sequences
 - ERCC annotations
 - a file containing a list of query SRA accessions
 
once you have these, the basic invocation is
```bash
nextflow run main.nf --ref-genome host.fa --genome-size $(wc -c host.fa | cut -d' ' -f1)
--ref-genome-gtf host.gtf --ercc-fa ERCC92.fa --ercc-gtf ERCC92.gtf --accession-list accession-list.txt -profile local,docker
```
or, by first creating a `params.yaml` file:
```yaml
ref-genome: host.fa
ref-genome-gtf: host.gtf
ercc-fa: ERCC92.fa
ercc-gtf: ERCC92.gtf
genome-size: 1396431182
accession-list: accession-list.txt
```

However, there are other parameters and configuration you will probably want.

First, note that `-profile local` causes the pipeline to use the `local` profile configuration in `nextflow.config`.  This is only appropriate for testing on a single user device.

If you are running on a cluster that runs slurm, this should be `-profile slurm` and that profile should be adapted to your cluster as appropriate (in particular, `maxMemory`, `maxCpus`, `queue` should be updated according to your partitions, and the fallback for the `scratch` temp directory should be customized)

If you are running on a cluster environment that uses a **different** workload manager/executor, you must copy the `slurm` profile to a new profile for your environment and, minimally, change the executor to the appropriate one according to https://www.nextflow.io/docs/latest/executor.html

Depending on your environment and workload, more configuration may be needed; see [[params.html]]


## On 'failed' Tasks
when running the pipeline you may get messages that some processes have failed.
```
[22/70ef62] process > prefetch (SRR7947912)          [ 99%] 185 of 186, failed: 1, retries: 1
[2d/76614c] process > fastq_dump (SRR11164691)       [ 89%] 164 of 184
[-        ] process > check_direct_fastqs            -
[2d/7cbd5b] process > filter_barcodes (SRR19077303)  [ 92%] 151 of 164
[20/299d1b] process > fastp (SRR16292403)            [ 96%] 148 of 153, failed: 2, retries: 2
[9d/aba73c] process > priceseqfilter (SRR12623494)   [ 89%] 131 of 146
[a4/1e902c] process > star_counts (SRR9608473)       [ 96%] 126 of 130
[cd/528a40] process > sort_bam (SRR7240614)          [ 99%] 125 of 126
[5d/ebdd65] process > htseq_count (SRR5687193)       [ 87%] 109 of 125
[9a/f711a0] process > hisat2 (SRR9608473)            [ 96%] 125 of 130
[16/5ca824] process > star (SRR6176732)              [100%] 125 of 125
[74/733caf] process > bowtie2 (SRR6268189)           [ 96%] 122 of 126, failed: 1, retries: 1
[bf/13a88e] process > process_bowtie2_sam (SRR120... [ 95%] 116 of 121
[29/c03fbb] process > bowtie2_filter_by_names (SR... [ 95%] 111 of 116
[d9/4eb1b9] process > dedup (SRR20001190)            [ 94%] 109 of 115, failed: 6, retries: 4
[71/ac0eb6] process > gsnap (SRR2968455)             [ 77%] 80 of 103
[01/d34f25] process > process_gsnap_sam (SRR5418456) [ 91%] 73 of 80
[24/fa8d08] process > gsnap_filter_by_names (SRR1... [ 95%] 70 of 73
[-        ] process > gsnap_skip                     -
[2c/30764f] process > stats_csv (SRR6933827)         [ 95%] 70 of 73
```
dedup failures are 'normal' when dealing with sequences from unknown sources -- if the average read
length is too short (by e.g. adaptor trimming or quality filtering) the sequence drops out due to
low information.

a small percentage of failures can also be expected from other steps due to "backoff retry" - when
the pipeline underestimates the memory required for an input, it will retry with increased memory
requests, up to a defined number of failures or maximum allowed memory request.  This behavior is
normal and usually benign (though it does make reading the progress summary more difficult).
SRA "prefetch" and "fastq-dump" are also known to be unusually prone to failure.  Throttling
concurrent downloads with `--parallel-downloads` can sometimes help alleviate this.

For other failures, especially large numbers of failures on *small* inputs, something may be wrong
worth investigating or tuning further.


## Troubleshooting
If STAR (or hisat2, or bowtie2, or gsnap) complain about missing input index files due to nextflow
dependency copy race conditions, you can try running 

`nextflow run modules/step.0.generate_indexes.nf`

to create the indexes ahead of time as a workaround.

Make sure to use the precomputed index [[parameters]] (`--hisat-ref-indexes`, `--star-ref-indexes-ercc`, `--star-ref-indexes`, etc.)