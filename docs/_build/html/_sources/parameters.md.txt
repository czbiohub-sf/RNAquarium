# RNAquarium pipeline parameters (Part I: Transcriptomic + Filtering)

## Reference genome options

Reference genome related files and options required for the workflow.

| Parameter | Description | Type | Default | Required |
|-----------|-------------|------|---------|----------|
| `genomeSize` | Approximate host genome size in bytes | `number` | — | Yes |
| `refGenome` | Host reference genome FASTA (bgzip compressed) | `string` | `Danio_rerio.GRCz11.dna_sm.primary_assembly.fa` | Yes |
| `refGenomeGtf` | Host reference GTF annotation | `string` | `Danio_rerio.GRCz11.108.gtf` | Yes |
| `erccFa` | ERCC spike-in FASTA sequences | `string` | `ERCC92.fa` | Yes |
| `erccGtf` | ERCC spike-in GTF annotation | `string` | `ERCC92.gtf` | Yes |
| `contamFa` | Contaminant FASTA for kallisto filtering | `string` | — | No |

:::{note}
The default genome filenames still reference GRCz11. The 75k production run used GRCz12tu (`GCF_049306965.1_GRCz12tu_genomic.fa.gz`) with `genome-size: 1448808683`, supplied via a params file that overrides these defaults.
:::

## Pre-built index options

Supply pre-built aligner indexes to skip index generation. If not provided, indexes are generated automatically from the reference genome.

| Parameter | Description | Type | Default |
|-----------|-------------|------|---------|
| `starRefIndexes` | STAR genome index directory (without ERCC) | `string` | — |
| `starRefIndexesErcc` | STAR genome index directory (with ERCC) | `string` | — |
| `hisatRefIndexes` | HISAT2 index directory | `string` | — |
| `bowtieRefIndexes` | Bowtie2 index directory | `string` | — |
| `gsnapRefIndexes` | GMAP/GSNAP index directory | `string` | — |
| `kbContamIndexes` | Kallisto contamination index file | `string` | — |

## Input/output options

| Parameter | Description | Type | Default | Required |
|-----------|-------------|------|---------|----------|
| `accessionList` | File with SRA run accessions (one per line, or CSV with `Run,size_MB` columns) | `string` | — | Yes (or `fastqPath`) |
| `fastqPath` | Glob pattern for pre-downloaded FASTQ files | `string` | — | No |
| `publishDir` | Output directory for results | `string` | `$PWD` | No |
| `tmp` | Scratch directory for temporary files | `string` | — | Recommended |
| `backupTmp` | Fallback scratch directory | `string` | — | No |

## Pipeline flow control

| Parameter | Description | Type | Default |
|-----------|-------------|------|---------|
| `skipHostCounts` | Skip the host gene counting branch (Step 3b) | `boolean` | `false` |
| `skipHisat` | Skip HISAT2 alignment step | `boolean` | `false` |
| `htseqCount` | Use HTSeq instead of featureCounts for gene counting | `boolean` | `false` |
| `outputCram` | Compress host alignments to CRAM format | `boolean` | `false` |
| `retainMixed` | Retain discordant pairs (one mate maps, one does not) as non-host | `boolean` | `true` |
| `kbRetainMixed` | Retain discordant pairs in kallisto filter | `boolean` | `false` |
| `hisatUseTranscript` | Include transcript/exon info in HISAT2 index | `boolean` | `true` |

## Resource and performance options

| Parameter | Description | Type | Default |
|-----------|-------------|------|---------|
| `parallelDownloads` | Maximum concurrent SRA downloads | `number` | `100` |
| `maxCpus` | Maximum CPUs per process | `number` | (profile-dependent) |
| `maxMemory` | Maximum memory per process | `string` | (profile-dependent) |
| `starThreadsSmall` | STAR threads for small files | `number` | `4` |
| `starThreadsLarge` | STAR threads for large files | `number` | `16` |
| `starSjdbOverhang` | STAR splice junction database overhang length | `number` | `100` |
| `starUseSharedMem` | Use STAR shared-memory genome loading | `boolean` | `false` |
| `seed` | Random seed for reproducibility | `number` | `32854` |

## Filtering options

| Parameter | Description | Type | Default |
|-----------|-------------|------|---------|
| `extraAdapters` | Path to additional adapter sequences FASTA | `string` | `$PWD/extra-adapters.fasta` |
| `sdFilterMates` | seq-detective: also filter the mate of a detected barcode file | `boolean` | `true` |
| `dedupPercentLen` | czid-dedup: percentage of read length threshold | `number` | `100` |
| `dedupMinLen` | czid-dedup: minimum read length for deduplication | `number` | `20` |

## Publishing options

Control which intermediate results are saved to the output directory.

| Parameter | Description | Type | Default |
|-----------|-------------|------|---------|
| `publishIntermediate` | Enable publishing of intermediate files | `boolean` | `false` |
| `publishFastqs` | Publish downloaded FASTQ files | `boolean` | `false` |
| `publishQCfiltered` | Publish fastp-filtered files | `boolean` | `false` |
| `publishReadcounts` | Publish host read count files | `boolean` | `true` |
| `publishKallisto` | Publish kallisto intermediate files | `boolean` | `false` |
| `publishHisat` | Publish HISAT2 intermediate files | `boolean` | `false` |
| `publishStar` | Publish STAR intermediate files | `boolean` | `false` |
| `publishBowtie` | Publish Bowtie2 intermediate files | `boolean` | `false` |
| `publishDedup` | Publish dedup intermediate files | `boolean` | `false` |
| `publishGsnap` | Publish GSNAP intermediate files | `boolean` | `true` |

## Advanced options

| Parameter | Description | Type | Default |
|-----------|-------------|------|---------|
| `cleanupIntermediate` | Automatically clean up intermediate files to save disk space | `boolean` | `true` |
| `nxfUnstageHack` | Enable workaround for Nextflow unstaging bug (see [Technical Notes](technical-notes.md)) | `boolean` | `false` |
| `backupScratchHack` | Use backup scratch directory for downloads | `boolean` | `false` |
