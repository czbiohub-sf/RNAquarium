The RNAquarium Nonhost Read Pipeline was developed by Chan Zuckerberg Biohub SF to identify viral
and other non-host sequences in publicly available Zebrafish datasets.

# Pipeline overview
## Step 0: Index generation
We use 4 mappers in the pipeline, all of which need genome indexes generated before mapping can
take place.

We generate two versions of the STAR indexes, one without ERCC spike-in controls as part of the
"reference genome" for use in the host counts pipeline (where we *do* want to map to host)

These steps take significant compute resources and time, but allow the mapping per-read to be fast.
Since we can start generating these in parallel during prefetch (below) and initial read filtering,
it doesn't translate to a strong bottleneck.

## Step 1: SRA prefetch and pre-processing
In this stage of the pipeline, we download the SRA accessions specified in `--accession-list`.

Because this download is slow and unreliable, we recommend that users with appropriate disk space
specify a `--fastq-path` **even when not using non-SRA database sequences** as this will store the
downloaded sequences and reuse them on subsequent pipeline runs rather than redownloading.

"fastq-dump" converts SRA files to FASTQ.

At this point, we have three patterns:
 1) single-ended reads, with one FASTQ file.
 2) paired reads, with two FASTQ files
 3) scRNA-seq reads, with one FASTQ file for barcode sequences and one FASTQ file for reads

The third case's two files are ambiguous with paired reads to the rest of pipeline, requiring
that we explicitly check whether a FASTQ looks like barcodes for scRNA-seq, and remove those.

At the same time as barcode characterization we collect various useful information about an
accession, such as read count, median read length, and file size in bytes.
These statistics are useful later for dynamic resource allocation.


## Step 2: Filtering and adapter trimming
In this step we do some basic processing of reads.
We use FASTP for adapter trimming and PriceSeqFilter for quality filtering.

## Step 3: HISAT2
Our first mapping step uses HISAT2