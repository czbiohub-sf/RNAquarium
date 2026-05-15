# Pipeline Overview (Part I: Transcriptomic + Filtering)

The RNAquarium transcriptomic and filtering pipeline was developed by Chan Zuckerberg Biohub to process publicly available zebrafish RNA-seq datasets. It downloads runs from the NCBI SRA, filters and trims reads, quantifies host gene expression, and passes reads through a series of host genome aligners to progressively remove host reads. The remaining unmapped reads are output as "possible non-host" for downstream metatranscriptomic analysis in Part II.

For the 75k zebrafish production run, this pipeline retained approximately 11 billion out of 1.64 trillion input reads (~0.7%) as unmapped non-host reads.

## Step 0: Index generation

The pipeline uses five alignment tools, each requiring its own set of genome indexes. If pre-built indexes are not provided via parameters, they are generated automatically at the start of the run:

| Aligner | Index type | Notes |
|---------|-----------|-------|
| HISAT2 | Splice-aware genome index | Optionally includes transcript/exon information (`--hisat-use-transcript`) |
| STAR | Genome index (without ERCC) | Used for host read counting |
| STAR + ERCC | Genome index (with ERCC spike-ins) | Used for the mapping/filtering step |
| Bowtie2 | Genome index | |
| GMAP/GSNAP | Hash-based genome index | |

Additionally, if a contaminant FASTA is provided (`--contam-fa`), a **kallisto** contamination index is built for the optional pseudoalignment-based host/contaminant filter.

Index generation is compute-intensive but runs in parallel with early pipeline steps (download, filtering) so it does not create a strong bottleneck. On subsequent runs, pre-built indexes can be supplied to skip this step entirely.

## Step 1: SRA download and pre-processing

*Nextflow processes: `download`, `check_direct_fastqs`, `filter_barcodes`*

Reads are downloaded from the NCBI SRA using `sra-tools` (prefetch + fasterq-dump). Because SRA downloads can be slow and unreliable, the pipeline supports throttling concurrent downloads (`--parallel-downloads`) and retries with backoff.

Alternatively, **local FASTQ files** can be provided via the `--fastq-path` parameter with a glob pattern. In this case the `check_direct_fastqs` process normalizes the Illumina-style filenames (`_R1_001.fastq.gz` / `_R2_001.fastq.gz`) to the internal naming convention expected by downstream steps, and the SRA download is skipped entirely.

During download, several pieces of metadata are collected per run: read count, median read length, file size in bytes, and single-end vs paired-end layout. These statistics drive dynamic resource allocation for all downstream steps.

For SRA inputs, **seq-detective** then analyzes each run to detect and filter technical sequences, particularly scRNA-seq barcodes. SRA datasets can contain three patterns: single-ended reads (one FASTQ), paired-end reads (two FASTQs), or scRNA-seq reads (one barcode FASTQ + one read FASTQ). The barcode files are ambiguous with paired-end layout, so seq-detective identifies and removes them. It can also optionally filter the corresponding mate file (`--sd-filter-mates`).

The `filter_barcodes` process then examines read lengths and discards barcode-only mates (median length < 32 bp with the other mate > 80 bp), converting such runs from paired-end to single-end.

Runs that produce empty or negligibly small files after filtering are dropped from further processing.

## Step 2: Adapter trimming and quality filtering

*Nextflow process: `fastp`*

**fastp** performs adapter trimming, quality filtering, and read-level QC. Custom adapter sequences can be provided via `--extra-adapters`. fastp also generates per-run HTML QC reports.

Reads that are too short or too low-quality after trimming are removed. Runs that fall below a minimum size threshold after this step are dropped.

## Step 2.5: Kallisto contamination filter (optional)

*Nextflow process: `kb_negative`*

If a kallisto contamination index is provided (`--kb-contam-indexes`), reads are pseudoaligned against it using **kallisto/bustools** via kb-python. Reads that align to known host or contaminant sequences are removed; unaligned reads continue to the mapping steps.

This optional step provides fast, lightweight host filtering before the more computationally expensive alignment steps. The `--kb-retain-mixed` parameter controls whether read pairs with one aligning mate are retained or discarded.

## Step 3: HISAT2

*Nextflow process: `hisat2`*

The first splice-aware alignment step uses **HISAT2** against the host genome. Reads that align to the host are removed; unmapped reads continue. For paired-end data, the `--retain-mixed` parameter controls whether discordant pairs (one mate maps, one does not) are kept or removed.

HISAT2 alignment statistics (aligned, multi-mapped, unmapped counts) are collected per run for the final summary.

## Step 3b: Host read counting (parallel branch)

*Nextflow processes: `star_counts`, `feature_count` (or `htseq_count`), `host_cram`*

In parallel with the non-host filtering path, a separate branch quantifies host gene expression. After fastp trimming, reads are aligned to the host genome with ERCC spike-ins using **STAR** in quantification mode. The resulting BAM is then processed by either:

- **featureCounts** (subread, default) — fast, multi-threaded gene counting
- **HTSeq** (optional, via `--htseq-count`) — stricter but slower counting

Per-run count vectors are collected into a combined `counts.csv` matrix. If `--output-cram` is enabled, STAR alignments are additionally compressed to CRAM format for archival.

Host counting can be skipped entirely with `--skip-host-counts`.

## Step 4: STAR

*Nextflow process: `star`*

Unmapped reads from HISAT2 are aligned to the host genome using **STAR**. Reads that align are removed; unmapped reads continue. STAR threading is dynamically adjusted based on file size (`--star-threads-small` / `--star-threads-large`) and available memory.

STAR can optionally use shared memory for the genome index (`--star-use-shared-mem`), which is useful when many concurrent STAR jobs run on the same node.

## Step 5: Bowtie2

*Nextflow process: `bowtie2`*

Unmapped reads from STAR are aligned using **Bowtie2** in its default end-to-end sensitive mode. This catches reads missed by the splice-aware aligners. Reads that align are removed.

Runs that produce empty output after this step are dropped before deduplication.

## Step 6: Deduplication

*Nextflow process: `dedup`*

**czid-dedup** removes PCR and optical duplicate reads. This step operates on the reads that survived all three alignment filters. Deduplication failures are common when average read length is very short (due to aggressive trimming), in which case the run is dropped with an `ignore` error strategy.

Parameters `--dedup-percent-len` and `--dedup-min-len` control the minimum length thresholds for deduplication.

## Step 7: GSNAP

*Nextflow process: `gsnap`*

The final alignment step uses **GSNAP** (from the GMAP suite) as a sensitive catch-all aligner. Reads that align to the host genome are removed. The remaining unmapped reads are the pipeline's final output — "possible non-host" reads written as compressed FASTQ files.

Output is organized into `unmapped_reads/Single/` and `unmapped_reads/Paired/` subdirectories by accession ID.

If GSNAP fails for a run (for reasons other than out-of-memory), the pipeline falls back to keeping the non-host reads from the previous step (post-dedup, pre-GSNAP).

## Summary statistics

*Nextflow processes: `stats_csv`, `cleanup_branched`*

At the end of the pipeline, per-run statistics from every step are collected into a CSV file (`stats-<timestamp>.csv`) and an overall summary (`host-filtering.summary.txt`). The CSV includes columns for read counts at each stage: starting reads, fastp output, kallisto output, HISAT2/STAR/Bowtie2 unmapped counts, dedup output, GSNAP unmapped counts, and final non-host read count.

## Output structure

```
<publish-dir>/
├── download/                    # Per-run download info and seq-detective results
├── host_counts/                 # Gene count matrices and alignments
│   ├── counts.csv               # Combined count matrix (all runs)
│   ├── Single/<id>/             # Per-run count files (single-end)
│   ├── Paired/<id>/             # Per-run count files (paired-end)
│   └── alignments/              # CRAM files (if --output-cram)
├── unmapped_reads/              # Final non-host reads (Part II input)
│   ├── Single/<id>/             # Per-run unmapped FASTQs (single-end)
│   └── Paired/<id>/             # Per-run unmapped FASTQs (paired-end)
├── reports/
│   └── stats-<timestamp>.csv   # Per-run statistics across all steps
├── host-filtering.summary.txt  # Aggregate summary
└── <index directories>/         # Generated aligner indexes (if not pre-built)
```
