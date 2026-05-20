# Technical Notes

## Nextflow Unstaging Hack
At the time when this pipeline was developed, Nextflow has a longstanding three-headed logic bug
in its handling of unstaging files from a scratch directory:
[nextflow-io/nextflow/issues/3711](https://github.com/nextflow-io/nextflow/issues/3711), [nextflow-io/nextflow/issues/4755](https://github.com/nextflow-io/nextflow/issues/4755), [nextflow-io/nextflow/issues/4718](https://github.com/nextflow-io/nextflow/issues/4718)

To discuss my solution we have to understand the nextflow `nxf_main` normal execution flow.

{lineno-start=1 emphasize-lines="4,6,7"}
```bash
nxf_main() {
	<...>
	set +e
	(set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
	pid=$!
	wait $pid || nxf_main_ret=$?
	nxf_unstage
	# 'afterScript' directive
}
```
the `nxf_launch` surrounded by other line noise runs the user process script, with various error
and logging options.

{lineno-start=5}
```bash
	pid=$!
	wait $pid || nxf_main_ret=$?
```
when it finishes, its exit code is saved to `nxf_main_ret`.
this variable is *the* source of the exit code.  it is used later in `on_exit()` to write to both
`.exitcode` and the nextflow wrapper's own return status.

*after* saving the exit code, `nxf_unstage` runs.  `nxf_unstage` **does not update nxf_main_ret**
if a problem occurs.  this is problem #1: nextflow does not propagate exit codes when unstaging.

problem #2 was identified by some Nextflow contributors in the issues listed above: 
in `nxf_unstage()`, `|| true` prevents stageout errors (e.g. from `nxf_fs_copy`) from propagating.
However, we actually don't mind this; it normally wouldn't be completely disastrous to continue
copying what does survive for the user to debug,
we want to survive to cleanup scratch space and run `afterScript`, et cetera.

problem #3 is that Nextflow does not provide a method to validate output after unstaging.
a 0-length `out.fastq` is just as valid as the well-formed result is just as valid as a file truncated
in the middle of the quality string.

### Manifestation
During testing, our lustre filesystem had a phase where some metadata accesses would go offsite and
fail to retrieve file information correctly.  When this (or other race conditions) occurred,
copying the output file from scratch to work directory would fail with
`cp: SRR12661676.trimmed.fastq: failed to get extents info: Invalid argument`.

In `nxf_unstage()`, each output file is being passed to `nxf_fs_copy`:

{emphasize-lines="7,9,10,11"}
```bash
nxf_unstage() {
	true
	cp .command.out	<workdir>/.command.out || true
	cp .command.err <workdir>/.command.err || true
	cp .command.trace <workdir>/.command.trace || true
	cp .command.env <workdir>/.command.env || true
	[[ ${nxf_main_ret:=0} != 0 ]] && return
	IFS=$'\n'
	for name in $(eval "ls -1d *.filtered.fastq.gz" | sort | uniq); do
		nxf_fs_copy "$name" <workdir> || true
	done
	unset IFS
}
```
{emphasize-lines="6"}
```bash
nxf_fs_copy() {
	local source=$1
	local target=$2
	local basedir=$(dirname $1)
	mkdir -p $target/$basedir
	cp -fRL $source $target/$basedir
}
```
problem 1: any error in `nxf_fs_copy` is not propagated to exitcode and hidden from the nextflow
runner.  errors are not handled at all other than (problem 2) to ignore them.

{start-line=7}
```[[ ${nxf_main_ret:=0} != 0 ]] && return```
notice the hint at problem 3: if the process script had a nonzero exit status, `nxf_unstage` avoids
unstaging whatever output files there might be, so that post-task processing doesn't see anything
that looks like valid output!

The result was that hundreds of runs would silently fail -- nextflow did not see them as errors
because one of the files in a mate pair was well-formed, or the filename existed but its contents 
could not be copied.  These runs continued through the pipeline and failed at unhelpful points,
already too late to send it back to the *true* failed step to retry.

### Solution
we need to replace nextflow's `nxf_fs_copy` with one that
 - checks for nonzero exit code from the `cp -fRL $source $target/$basedir`
 - waits before trying again, for networked filesystem race conditions
 - if still unable to unstage, *deletes outputs so that nextflow detects the failure properly*

```bash
nxf_fs_copy() {
	local source=$1
	local target=$2
	local basedir=$(dirname $1)
	mkdir -p $target/$basedir
	cp -fRL $source $target/$basedir
	nxf_fs_copy_ret=$?
	if [[ $nxf_fs_copy_ret > 0 ]]
	then
		sleep 10  # second chance
		cp -fRL $source $target/$basedir
		nxf_fs_copy_ret=$?
		if [[ $nxf_fs_copy_ret > 0 ]]
		then
			ls -alh
			echo "NXF_UNSTAGE copy $source failed ($nxf_fs_copy_ret)" | tee .command.err
			local temp=$(basename $source)
			rm -fr $target/$basedir/${temp//_[12].fastq/_*.fastq} 
		fi
	fi
}
```


## Running with local FASTQs

Instead of downloading from SRA, you can provide local FASTQ files using the `fastq-path` parameter with a glob pattern:

```yaml
fastq-path: /path/to/fastqs/*_R{1,2}_001.fastq.gz
```

The expected filename convention is Illumina-style: `{sampleID}_R1_001.fastq.gz` and `{sampleID}_R2_001.fastq.gz` for paired-end data. Nextflow's `fromFilePairs` extracts the sample ID as everything before `_R1_001` or `_R2_001`. For example, `chikv_1_S9_R1_001.fastq.gz` produces sample ID `chikv_1_S9`. These IDs appear in all output filenames and directory names.

Single-end reads and `.fq.gz` variants are also supported. A flat directory of FASTQs works fine (no subfolders needed).

When using `fastq-path`, omit `accession-list` from your params file. The pipeline will use the `check_direct_fastqs` process instead of the SRA download process.

Here is a complete example `params.yaml` for local FASTQs with the GRCz12tu zebrafish reference:
```yaml
# Input
fastq-path: /path/to/fastqs/*_R{1,2}_001.fastq.gz

# Reference genome (GRCz12tu)
ref-genome:     /path/to/GCF_049306965.1_GRCz12tu_genomic.fa.gz
ref-genome-gtf: /path/to/GCF_049306965.1_GRCz12tu_genomic.gtf
genome-size: 1448808683

# ERCC spike-ins
ercc-fa:  /path/to/ERCC92.fa
ercc-gtf: /path/to/ERCC92.gtf

# Output
publish-dir: output
tmp: /tmp/

# Pipeline behavior
output-cram: false
cleanup-intermediate: true
retain-mixed: true
```


## Reference index directory names

The default output directory names for mapper indexes are:

 - `hisat2_{Genus_species}_genome`
 - `star_{Genus_species}_indexes.ERCC`
 - `star_{Genus_species}_indexes`
 - `bowtie2_{Genus_species}_index`
 - `gmap_{Genus_species}_genome`

Where `{Genus_species}` comes from the first part of the reference genome filename before a dot (`.`), e.g., `Danio_rerio` from `Danio_rerio.GRCz11.dna_sm.primary_assembly.fa`, or `GCF_049306965` from `GCF_049306965.1_GRCz12tu_genomic.fa.gz`.

These names can be provided to the reference index parameters (`--hisat-ref-indexes`, `--star-ref-indexes-ercc`, `--star-ref-indexes`, `--bowtie-ref-indexes`, `--gsnap-ref-indexes`) for subsequent runs of the pipeline, saving significant time on re-runs.


## Utility Scripts

The following utility scripts are provided in `src/nonhost/` for pipeline management and debugging:

**`delete-step.sh`** — Removes all work directories for a specific Nextflow process, allowing that step to be re-run from scratch on the next `-resume`. Usage: pass the process name as an argument (e.g., `bash delete-step.sh gsnap`). This is useful when a step needs to be re-executed with different parameters without re-running the entire pipeline.

**`collect-reruns.sh`** — Identifies SRA accessions that failed to complete the pipeline by comparing the input accession list against completed output directories. It automatically excludes expected dropouts (dedup failures due to short reads). Outputs the list of accessions that need re-running, either as a plain text list or as a filtered RunInfo CSV. Usage: `bash collect-reruns.sh input_accessions.txt nonhost_reads_dir counts_dir .nextflow.log`

**`util/debug_run.sh`** — Diagnoses failures from `.nextflow.log` by extracting all accessions that terminated with errors, then for each one printing the process name, SLURM job ID, working directory, `seff` resource usage summary, and the task's `.command.log`. This produces a per-accession log file useful for diagnosing whether failures were due to OOM, preemption, or other causes. Usage: `bash util/debug_run.sh .nextflow.log`

**`util/dropouts.sh`** — Generates a summary of pipeline dropouts and failures by category: download failures, corrupt files, dedup failures, gsnap failures (skipped runs), preemptions, and other errors. Outputs a single CSV line with counts. Usage: `bash util/dropouts.sh nonhost_reads_dir .nextflow.log`. Set `VERBOSE=1` for per-failure details.

**Submission script templates** (`nextflow-submit-example.sh`, `nextflow-submit-local-test.sh`) — SLURM submission script templates showing recommended `#SBATCH` directives, module loading, and Nextflow invocation patterns.


## Dependencies

There are multiple ways to obtain the dependencies for the RNAquarium preprocessing pipeline:

 - **Conda/mamba packages** (preferred for ease of deployment) — conda packages may not have optimizations for the user's hardware. Some tools are not available as conda packages, and `setup-minimal.sh` is provided for these.
 - **Local installation** — if it works, it works. May be unsafe for heterogeneous computing clusters, as no attempt is made to identify different machine types.
 - **A singular "RNAquarium" container** with all tools — may not have all optimizations for the user's hardware. May have limited network bandwidth, which slows down the prefetch step.
 - **Individual biocontainer packages** — not all packages used are available, and `setup-minimal.sh` is provided for these. May have limited network bandwidth, which slows down the prefetch step.

Some dependencies present additional problems which require careful testing when using a container. Newer versions of gmap/gsnap are much faster at mapping, but have a bad container, different sensitivity, fail to detect all the input reads, or crash on certain input sequences.

```
gsnap_version,         gsnap_time,
after_dedup,total_read,aligned,multialign,aligned_unique,unaligned,mixed
2021.08.25--pl5321h...,     99.0s,
      17507,     17507,   1221,         0,          1221,    16286,    0
2021.12.17(self-compiled),  27.4s,
      17507,     17507,   1150,         0,          1150,    16357,    0
! 2023.03.24--h67092d7_0 doesn't work (bad hardcoded perl path in container)
! 2023.03.24--h9d449c0_1 doesn't work (bad hardcoded perl path in container)
! 2023.04.28--h9d449c0_1 doesn't work (bad hardcoded perl path in container)
! 2023.07.20--h9d449c0_1 doesn't work (bad hardcoded perl path in container)
2023.10.10--h9d449c0_1,     16.1s,
      17507,      2399,    158,         0,           158,     2241,    0
2023-12-01(self-compiled),   6.1s,
      17507,     17507,    275,         0,           275,    17232,    0
! 2024-02-22(self-compiled) inconsistent (segfaults on some inputs)
```
Other dependencies may or may not be safe to up- or downgrade.
