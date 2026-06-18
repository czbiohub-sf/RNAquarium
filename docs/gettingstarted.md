# Getting Started

## Dependencies

The RNAquarium pipeline requires several tools to be available in `$PATH`.

They can either be installed separately, or, more typically, provided by containers or conda environments. Use `-profile docker`, `-profile singularity`, or `-profile conda`/`-profile mamba` as appropriate for your environment.

There are a small number of programs that are not available as containers or conda packages. These can be installed with:
```bash
cd src/nonhost
bash setup-minimal.sh
```
This installs `fastq-lengths`, `fastq-namefilter`, `fastq-numfilter`, `PriceSeqFilter`, `czid-dedup`, and `gsnap`/`gmap` into `src/nonhost/bin/`. Make sure this directory is in your `$PATH` when running the pipeline (the submission script templates handle this automatically).

### `seq-detective` (required only for SRA downloads)

`seq-detective` is required **only** when running from SRA accessions. It is invoked exclusively by the `download` step, where it subsamples the downloaded reads, maps them against the host index, and classifies each mate as biological or technical (e.g. a single-cell cell-barcode read, a poly-A/poly-T tail, or a low-information mate), removing technical mate files before downstream processing. If you run **only from local FASTQ files** you do not need it: the local path uses a lighter, read-length-based check (`fastq-lengths`, installed above) instead. No other step uses `seq-detective`.

`seq-detective` is **not** installed by `setup-minimal.sh`. Install it with:
```bash
cd src/nonhost
nextflow run install-deps.nf
```
This clones [seq-tech-detective](https://github.com/czbiohub-sf/seq-tech-detective) and builds `seq-detective` into `src/nonhost/bin/`. The build uses `meson` and a C/C++ compiler with OpenMP, plus `mamba`/`conda` for runtime dependencies — these are provided by the repo's `environment.yml`, so ensure `mamba` (or `conda`) is available (e.g. `module load mamba`).

**Install layout:** the build produces three directories that must be kept together under `src/nonhost/` — `bin/` (the `seq-detective` launcher), `libexec/seq-detective/` (the subcommands), and `share/seq-detective/` (barcode onlists). The launcher locates `libexec` and `share` relative to its own path, so keep all three as siblings; copying only the `bin/` binary will cause runtime "subcommand not found" errors.

> **Note:** on clusters where compute nodes have no outbound network or git access, the automated clone in `install-deps.nf` may not run from a batch job. In that case, build `seq-detective` once on a login node and copy the three directories into `src/nonhost/`.


## Running the pipeline

The main preprocessing pipeline is the default entry point of `main.nf` inside `src/nonhost/`.

### Minimum required inputs

 - A host genome (bgzip-compressed `.fa.gz` with `.gzi` and `.fai` indexes)
 - Host genome annotations (`.gtf`)
 - An estimate of the host genome size in bytes
 - ERCC spike-in sequences and annotations
 - Either an SRA accession list **or** a glob pattern pointing to local FASTQ files

### Running from SRA accessions

```bash
nextflow run main.nf --ref-genome host.fa.gz --genome-size 1448808683 \
  --ref-genome-gtf host.gtf --ercc-fa ERCC92.fa --ercc-gtf ERCC92.gtf \
  --accession-list accession-list.txt -profile slurm,mamba
```

Or by creating a `params.yaml` file:
```yaml
ref-genome: host.fa.gz
ref-genome-gtf: host.gtf
ercc-fa: ERCC92.fa
ercc-gtf: ERCC92.gtf
genome-size: 1448808683
accession-list: accession-list.txt
```

## Preparing for Part II (metatranscriptomics)

If you plan to run the metatranscriptome pipeline on these samples afterward, you will need a `bioproject_mapping.json` that groups sample IDs into assembly groups. The group IDs must start with "PRJ" followed by two uppercase letters and digits (e.g., `PRJWILD00001`), since the Part II code validates this format with a regex. See [Inputs and Databases](inputs-and-databases.md) for details.


## HPC profiles and customization

The pipeline includes several profiles defined in `nextflow.config`:

 - **`local`** — for testing on a single machine; not suitable for production runs.
 - **`slurm`** — generic SLURM profile; a good starting point for most HPC clusters.
 - **`bruno`** — CZ Biohub-specific profile with site-specific queue names and options. External users should use `slurm` instead or customize `bruno` for their site.
 - **`docker`**, **`singularity`** — container runtimes.
 - **`conda`**, **`mamba`** — conda-based dependency management.

Typical usage combines an executor profile with a dependency profile, e.g., `-profile slurm,mamba`.

If you need to customize SLURM options (partition names, QoS, temp disk reservations), you have two options:

1. **Edit `nextflow.config` directly** — simplest approach, but be aware that Nextflow's `-resume` caches task scripts including SLURM parameters. If you change config after a run, you must purge both `.nextflow/` and the work directory for changes to take effect.
2. **Use an override config file** — create an `override.config` and pass it with `-c override.config`. Note that with Nextflow 23.10.x, `-c` overrides may not take effect when combined with `-resume` on cached tasks.


## Submission script

The recommended way to launch the pipeline on an HPC is via a SLURM submission script. The key points:

 - **Run `sbatch` from the repository root directory.** The script uses `SLURM_SUBMIT_DIR` to locate the repo and handles `cd src/nonhost` internally. This is necessary because Nextflow interprets `main.nf` as a local file only when run from the directory containing it; otherwise it may try to fetch it as a GitHub repository URL.
 - **The script adds `src/nonhost/bin/` to `$PATH`** so locally compiled tools (from `setup-minimal.sh`) are available.

See the submission script templates in `src/nonhost/` for examples.


## On 'failed' Tasks

When running the pipeline you may see messages that some processes have failed.
```
[22/70ef62] process > download (SRR7947912)          [ 99%] 185 of 186, failed: 1, retries: 1
[-        ] process > check_direct_fastqs            -
[2d/7cbd5b] process > filter_barcodes (SRR19077303)  [ 92%] 151 of 164
[20/299d1b] process > fastp (SRR16292403)            [ 96%] 148 of 153, failed: 2, retries: 2
```

Dedup failures are normal when dealing with sequences from unknown sources — if the average read length is too short (from adapter trimming or quality filtering), the sequence drops out due to low information content.

A small percentage of failures can also be expected from other steps due to backoff retry: when the pipeline underestimates the memory required for an input, it retries with increased memory requests, up to a defined number of failures or maximum allowed memory. This behavior is normal and usually benign.

On preemptible partitions (e.g., `preempted`), exit status 143 (SIGTERM) indicates the job was preempted by a higher-priority job. The pipeline retries these automatically.

SRA "prefetch" and "fasterq-dump" are also known to be unusually prone to failure. Throttling concurrent downloads with `--parallel-downloads` can sometimes help.

For other failures, especially large numbers of failures on small inputs, something may be worth investigating or tuning further. The utility script `util/dropouts.sh` can help identify which samples dropped out and at which step.


## Troubleshooting

If STAR (or HISAT2, Bowtie2, or GSNAP) complains about missing input index files due to Nextflow dependency copy race conditions, you can try running

```bash
nextflow run modules/local/step.0.generate_indexes.nf
```

to create the indexes ahead of time as a workaround. Make sure to use the precomputed index [parameters](parameters.md) (`--hisat-ref-indexes`, `--star-ref-indexes-ercc`, `--star-ref-indexes`, etc.) on subsequent runs.

**Nextflow caching pitfall:** If you change `nextflow.config` or `override.config` after a previous run, Nextflow's `-resume` may continue using cached task scripts with the old settings. To force regeneration, delete both the `.nextflow/` cache directory and the work directory, then re-run without `-resume`.
