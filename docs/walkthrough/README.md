# RNAquarium Walkthrough — First-Time Test Run on an HPC

A step-by-step guide to test-running RNAquarium on a SLURM cluster with a tiny,
fully reproducible two-sample zebrafish example. It starts from an empty directory
and a fresh conda sandbox and assumes little prior `git`, `conda`, or Nextflow.

**Covers Part I (Transcriptomic + Filtering / "nonhost") and the Nextflow phase of
Part II (Metatranscriptomics), each validated end to end on a non-Biohub cluster.**
Part II's downstream (combine, curation, Salmon quantification) runs as separate
SLURM/R scripts and is out of scope here.

Example dataset: **BioProject PRJEB28062**, runs **ERR2723539** (~11 MB) and
**ERR2723409** (~38 MB) — small enough to run fast, real enough to reproduce.

---

## What's in this folder

In the repository these files live under **`docs/walkthrough/`**. The table is in the
order you use them.

**Everything here is a read-only template.** Copy what you need into your own run area
(or into the repo, where the table says so) and edit *the copy* — never the original.
A later `git pull` then can't clobber your edits, and your outputs never get committed.

| File | Action | Used at | Where it goes, and what you change |
| --- | --- | --- | --- |
| `README.md` | read | start | you are here |
| `zf-example.csv` | copy as-is | Part I · Step 4 | your run area — nothing to change |
| **`params.zf-example.yaml`** | **EDIT** | Part I · Step 5 | your run area — every `/path/to/…`, plus `genome-size` |
| **`nextflow-submit-part1-nonhost.sh`** | **EDIT** | Part I · Steps 6–7 | your run area — `--partition`, `PROJECT_DIR`, the conda Case A/B/C lines, `--genome-size` |
| `nextflow.config.walkthrough` | copy as-is | Part I · Step 6 | **into the repo**, over `src/nonhost/nextflow.config` |
| `stats.walkthrough.nf` | copy as-is | Part I · Step 6 | **into the repo**, over `src/nonhost/modules/local/stats.nf` |
| `step.7.gsnap.walkthrough.nf` | copy as-is | Part I · Step 6 | **into the repo**, over `src/nonhost/modules/local/step.7.gsnap.nf` |
| `zf-example-bioproject-mapping.json` | copy as-is | Part II · Step 3 | your run area — nothing to change |
| **`build-toy-databases.sh`** | **EDIT** | Part II · Step 4 | your run area — `DB_DIR`, optionally the seeded virus accessions. Toy-DB route only |
| `modules.walkthrough.config` | copy as-is | Part II · Step 4 | **into the repo**, over `src/metatranscriptome/conf/modules.config` |
| **`run-part2-metatranscriptome.sh`** | **EDIT** | Part II · Step 5 | your run area — `--partition`, the three literal `#SBATCH` paths, `PROJECT_DIR`, `DB_DIR` |
| `part2-test-databases.md` | read (optional) | Part II · Step 4 | background for the two toy-DB rows above |

The four files copied **into the repo** are the "swaps". They are `git`-revertible
(`git checkout -- <path>`), and each step below shows how to back them up first.

---

## The mental model

Three pieces cooperate: **Nextflow** (orchestrator — reads `main.nf`, builds the step
graph, submits each step as a job), **SLURM** (scheduler — runs each step on a compute
node), and **conda/mamba or Singularity** (supplies each step's tools automatically via
`-profile mamba` for Part I / `-profile singularity` for Part II).

**Part I produces** (per dataset) a host gene-counts table and a set of non-host reads.
**Part II** consumes those non-host reads.

---

## Setup — Prerequisites

Run each; you want a version number, not "command not found":

```bash
sbatch --version      # SLURM scheduler
java -version         # Nextflow needs Java 11+ (17 is fine)
git --version
```

**Networking:** the **login node** usually has internet but compute nodes may not — do
downloads and one-time builds on the login node.

### Getting conda working (do this once — it sets commands you'll use later)

How you make conda available differs by cluster. Determine which case you're in — the
"SITE SETUP" block in the submission scripts must match it:

```bash
conda --version                                       # if this works, likely Case C
module avail 2>&1 | grep -iE 'conda|mamba|anaconda'   # conda-ish modules, if any
```

- **Case A — conda from a module** (common): `conda --version` fails until you
  `module load <name>` (e.g. `module load anaconda`).
- **Case B — personal Miniforge/Miniconda**: no module; your line is
  `source ~/miniforge3/etc/profile.d/conda.sh` (adjust path).
- **Case C — conda always on PATH**: no setup line needed.

A SLURM batch job runs a **fresh shell that doesn't inherit** your interactive conda,
so the scripts re-establish it using your Case A/B/C line. Remember your letter.

---

## Setup — Conda environments (two of them)

RNAquarium's two parts use **different Nextflow versions**. At Biohub these come from
`module load nextflow/...`; external users won't have those modules, so create one
small conda env per version:

```bash
# Part I: Nextflow 23.10.1 + mamba 1.5.x (see Known Issues for why pinned) + samtools/htslib (genome prep)
conda create -n rnaq-part1 -c conda-forge -c bioconda nextflow=23.10.1 "mamba=1.5.*" samtools htslib

# Part II: Nextflow only — its pipeline tools run inside Singularity containers (samtools as a convenience)
conda create -n rnaq-part2 -c conda-forge -c bioconda nextflow=24.10.5 samtools
```

Confirm the Part I env got the real mamba (must be **1.5.x**, not `0.1.2`):

```bash
conda activate rnaq-part1
nextflow -version          # 23.10.1
mamba --version            # mamba 1.5.x
```

**Optional third env — Part II toy databases** (only if you do the Part II test-DB
step; see `part2-test-databases.md`). Needs BLAST + DIAMOND + entrez-direct + the R
`taxonomizr` package, all from conda:

```bash
conda create -n rnaq-dbbuild -c conda-forge -c bioconda blast diamond entrez-direct r-taxonomizr
conda activate rnaq-dbbuild
```

(If conda's R misbehaves on your cluster, you can use your own R with `taxonomizr`
instead — see Known Issues.)

| Env | Purpose | Contents |
| --- | ------- | -------- |
| `rnaq-part1` | run Part I | `nextflow=23.10.1`, `mamba=1.5.*`, `samtools`, `htslib` |
| `rnaq-part2` | run Part II | `nextflow=24.10.5`, `samtools` (tools come from containers) |
| `rnaq-dbbuild` | build Part II toy DBs (optional) | `blast`, `diamond`, `entrez-direct` + R/`taxonomizr` |

> The metatranscriptome README lists SPAdes/BLAST/DIAMOND/BBTools/Salmon/etc. as
> "required in PATH" — that's for running **without** containers. With
> `-profile singularity`, Nextflow supplies all of them; you don't install them.

**Which environment to activate at each step** (each step below also states this at the
top, so you never have to remember — but here's the whole map):

| Step | `conda activate …` |
| ---- | ------------------ |
| Part I · Steps 1–8 (setup-minimal, install-deps, refs, edits, the `sbatch` run) | **`rnaq-part1`** |
| Part II · Step 4 — building the toy databases | **`rnaq-dbbuild`** |
| Part II · Step 4 — `apptainer build` the container | **`conda deactivate`** (clean shell) |
| Part II · Steps 5–6 (edit + the `sbatch` run) | **`rnaq-part2`** |

Every environment is active only in *your login shell*; the `sbatch` submission scripts
re-activate the right one *inside the batch job* themselves, so the activation you do
interactively is for the hand-run commands (and to keep things unambiguous).

---

## Setup — Clone and orient

```bash
git clone https://github.com/czbiohub-sf/RNAquarium
```

Part I lives in `src/nonhost/` (its `main.nf`, `nextflow.config`, `setup-minimal.sh`,
`install-deps.nf`); Part II in `src/metatranscriptome/`. You run each pipeline's
`main.nf` from *its own* directory.

---

## Where these files go (and how to run them)

Keep the repo, your run area, and scratch separate. Set two shell variables to your
own paths, then copy the templates across:

```bash
REPO=/path/to/RNAquarium          # your git clone
RUN=/path/to/your/run-area        # where you edit + submit from (NOT inside the repo)

# Part I inputs + submission script:
cp "$REPO"/docs/walkthrough/zf-example.csv                    "$RUN"/
cp "$REPO"/docs/walkthrough/params.zf-example.yaml            "$RUN"/
cp "$REPO"/docs/walkthrough/nextflow-submit-part1-nonhost.sh  "$RUN"/

# Part II (later): mapping + optional toy-DB builder + Part II script
cp "$REPO"/docs/walkthrough/zf-example-bioproject-mapping.json "$RUN"/
cp "$REPO"/docs/walkthrough/run-part2-metatranscriptome.sh     "$RUN"/
cp "$REPO"/docs/walkthrough/build-toy-databases.sh             "$RUN"/   # toy-DB run only
```

The four swaps go *into the repo* instead, at their steps below (Part I Step 6, Part
II Step 4) — see the table above for their destinations.

---

> **Before the Part I steps: `conda activate rnaq-part1`, and keep it active.** Steps 2
> and 5–7 need `nextflow` (and, in the batch job, `mamba`) on your `PATH`, which this env
> provides. Step 1 (compiling the bundled tools) doesn't strictly need it, but activating
> now and staying in the env through Part I is simplest.

## Part I · Step 1 — Install the pipeline's local dependencies  (~10 min)

From `src/nonhost/`, build the small bundled tools:

```bash
cd "$REPO"/src/nonhost
bash setup-minimal.sh 2>&1 | tee setup-minimal.log
echo "exit code: ${PIPESTATUS[0]}"                       # want 0
```

Confirm the seven tools landed and are executable — **this is the real check**:

```bash
# Anchor on an absolute path -- "bin/" is relative, so running this from the wrong
# directory reports seven false MISSes rather than a missing directory.
BIN="$REPO/src/nonhost/bin"          # or type the full path to your clone
if [ ! -d "$BIN" ]; then
  echo "WRONG PATH: $BIN does not exist — check \$REPO, or that setup-minimal.sh ran"
else
  for t in fastq-lengths fastq-namefilter fastq-numfilter PriceSeqFilter czid-dedup gsnap gmap; do
    if [ -x "$BIN/$t" ]; then echo "OK   $t"; else echo "MISS $t"; fi
  done
fi
```

Exit code `0` + seven `OK`s means you're set. (If `fastq-numfilter` is missing, see
Known Issues.)

> **If you get seven `MISS`es, check where you are before suspecting the build.**
> These are files on disk, not commands on `PATH` — the check consults neither `PATH`
> nor your conda env, so an active `rnaq-part1` makes no difference. All seven missing
> at once almost always means the path is wrong, not that seven builds failed.

---

## Part I · Step 2 — Build `seq-detective` (SRA route only)  (~5 min)

```bash
# from src/nonhost/
nextflow run install-deps.nf
```

(If it fails with `prefix already exists: .../seq-detective`, see Known Issues.)

---

## Part I · Step 3 — Reference inputs

Download the host genome (GRCz12tu) from NCBI Datasets — the
[GCF_049306965.1 page](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_049306965.1/).
On the download dialog, **tick both "Genome sequences (FASTA)" AND "Annotation features
→ Genomic coding sequences" so you get the GTF** — Part I needs the annotation for host
gene counting, and a FASTA-only download will fail later. Make sure the annotation format
is **GTF** (not GFF). Save the zip into a `refs/` folder, then:

```bash
cd /path/to/your/refs
unzip ncbi_dataset.zip
ls ncbi_dataset/data/GCF_049306965.1/
#   GCF_049306965.1_GRCz12tu_genomic.fna   <- genome FASTA (uncompressed)
#   genomic.gtf                            <- annotation (GTF — this must be present)
```

> If you don't see `genomic.gtf` (only the `.fna`), the annotation wasn't included in the
> download — go back and re-download with the annotation/GTF box ticked.

> ⚠️ The NCBI zip contains its own `README.md` — if you unzip where a `README.md`
> exists, answer **`n`** to the overwrite prompt. Keeping refs in a separate folder
> avoids this.
> If the zip has both `GCA_...` and `GCF_...`, use **GCF** (it ships the GTF).

Only the genome FASTA must be bgzip-compressed + faidx-indexed; **it must end in
`.fa.gz` (not `.fna.gz`)** — so rename during copy:

```bash
cp ncbi_dataset/data/GCF_049306965.1/GCF_049306965.1_GRCz12tu_genomic.fna  GCF_049306965.1_GRCz12tu_genomic.fa
cp ncbi_dataset/data/GCF_049306965.1/genomic.gtf  GCF_049306965.1_GRCz12tu_genomic.gtf
bgzip GCF_049306965.1_GRCz12tu_genomic.fa
samtools faidx GCF_049306965.1_GRCz12tu_genomic.fa.gz
```

**Genome size** (bp) — use your file's real value (feeds resource scaling only):

```bash
awk '{s+=$2} END{print s}' GCF_049306965.1_GRCz12tu_genomic.fa.gz.fai
```

**ERCC spike-ins:**

```bash
cd /path/to/your/refs
wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip
unzip ERCC92.zip                 # -> ERCC92.fa  and  ERCC92.gtf
```

---

## Part I · Step 4 — The input list

`zf-example.csv` is a RunInfo-style CSV with a header and one accession per line
(`Run,size_MB`); the `size_MB` column drives resource scaling. It's ready to use as-is
for the example.

---

## Part I · Step 5 — Edit `params.zf-example.yaml`

Open **your copy** and replace every `/path/to/...`. Parameter names are **kebab-case**
(`ref-genome`, `ref-genome-gtf`, `genome-size`, `publish-dir`; `publishQCfiltered` is
the one camelCase exception). Set `genome-size:` to your `awk` value and `publish-dir:`
to your output dir. The pre-built-index block stays commented (indexes auto-generate on
first run); `contam-fa` is optional.

---

## Part I · Step 6 — Edit the submission script and swap in the test config

**a) Edit `nextflow-submit-part1-nonhost.sh`** (your copy):

- `PROJECT_DIR` — your run area (`SCRATCH_DIR`, `PARAMS_FILE`, `WORKDIR`, `LOGDIR`,
  `REPORTDIR` derive from it; leave those).
- `--partition` in the `#SBATCH` header — your cluster's partition name.
- `RESUME` — leave `"no"`. A fresh run refuses to write on top of an existing
  `unmapped_reads/`; set `"yes"` only to continue an interrupted run using the same
  work dir.
- `EXTRA_CONFIG` — leave empty. It exists to load a diagnostic config overlay.
- **SITE SETUP block** — two steps: (2a) make conda available with your Case A/B/C line
  (`module load anaconda` / `source ~/…/conda.sh` / nothing); (2b) the
  `source "$(conda info --base)/…/conda.sh"` + `conda activate rnaq-part1` lines always
  run. Both A and B still need 2b. (On many clusters A + 2b are both used — normal.)
- **`--genome-size`** in the run block — edit to your genome's `awk` value.

**b) Swap in the three test-scale Part I files** (all `git`-revertible):

```bash
cd "$REPO"
# back up the originals first (independent of git, so they survive a branch switch)
for f in src/nonhost/nextflow.config \
         src/nonhost/modules/local/stats.nf \
         src/nonhost/modules/local/step.7.gsnap.nf ; do
    [[ -e "$f.bak" ]] || cp -p "$f" "$f.bak"
done

cp docs/walkthrough/nextflow.config.walkthrough   src/nonhost/nextflow.config
cp docs/walkthrough/stats.walkthrough.nf          src/nonhost/modules/local/stats.nf
cp docs/walkthrough/step.7.gsnap.walkthrough.nf   src/nonhost/modules/local/step.7.gsnap.nf
# revert later:
#   git checkout -- src/nonhost/nextflow.config \
#                   src/nonhost/modules/local/stats.nf \
#                   src/nonhost/modules/local/step.7.gsnap.nf
```

What each one does:

- **`nextflow.config`** — strips Biohub-specific `--qos`/`--tmp`, sets conda channels,
  guards gsnap's `cpus` directive against a null `ext.largeThreshold`, and sets
  `gsnap`/`host_cram` `publishDir` to `mode: 'copy'` instead of `'move'`. `move`
  removes a process's outputs from the work dir, which makes `-resume` unable to
  cache those steps and leaves nothing recoverable if a publish fails.
- **`stats.nf`** — prevents an end-of-run summary crash on tiny inputs. Changes no
  counts or math.
- **`step.7.gsnap.nf`** — makes gsnap failures *detectable*. Upstream tests
  `if [[ $? > 0 ]]` after a `set -e`, which reads `set -e`'s status rather than
  gsnap's **and** string-compares rather than numerically — so it is never true and
  gsnap's "pass the reads through unfiltered" fallback can never run. See
  Known Issues → "Part I reports success but produces no non-host reads".

All three are **test-only** — revert for real large-data runs.

---

## Part I · Step 7 — Run via `sbatch`  (~25 min first run; cached after)

Submit **from the repo root** (the script uses `SLURM_SUBMIT_DIR` to find the repo):

```bash
conda activate rnaq-part1             # the Part I env (harmless if already active)
cd "$REPO"
mkdir -p nonhost_slurm.out            # SLURM needs this before submit
sbatch "$RUN"/nextflow-submit-part1-nonhost.sh
squeue --me
```

Watch progress (the log is timestamped per run, so pick the newest):

```bash
tail -f "$(ls -t "$RUN"/logs/nextflow_part1_zf-example-*.log | head -1)"
```

The first run builds aligner indexes (~25 min). Later runs can reuse them, but only
if you set `RESUME="yes"` in the script — it is `"no"` by default so that a re-test
is genuinely clean and never silently inherits a previous run's cache.

Each run also writes its own `reports/part1-<timestamp>/` directory in your run area
containing `trace-*.txt`, which records the status and exit code of every task. That
trace is the first thing to read if anything looks wrong.

---

## Part I · Step 8 — Find your outputs, and verify them

Under your `publish-dir`, the key Part I deliverable is
`unmapped_reads/Paired/<accession>/` (the non-host reads Part II consumes), plus host
count tables and a `reports/` summary CSV.

**Verify before moving on — the directories existing is not enough.** Part I can print
`Execution complete -- Goodbye` and exit 0 while having produced no reads at all: the
gsnap step's `errorStrategy` resolves to `'ignore'`, so a failed gsnap is dropped from
the channel without failing the workflow. The submission script checks this for you and
exits non-zero with diagnostics if anything is missing, but to confirm by hand:

```bash
PUB=$(sed -n 's/^publish-dir:[[:space:]]*//p' "$RUN"/params.zf-example.yaml)
find "$PUB/unmapped_reads" -name '*.fastq.gz' -size +132c -printf '%10s  %p\n' | sort -k2
```

You should see **two `.fastq.gz` files per sample**, each at least tens of KB — plus
`gsnap_out.sam` and a `gsnap.stats.txt` of a few hundred bytes. For the two-sample
example that is four FASTQs in total.

Two things to look at in the filenames:

- **`...gsnapFiltered.fastq.gz`** — normal. gsnap ran and the host-mapped reads were
  filtered out.
- **`...gsnapSkipped.fastq.gz`** — gsnap **failed**, and the pipeline passed the
  upstream reads through unfiltered. Part II will still run, but the reads have not had
  the final host-filtering pass. Check the gsnap step's `.command.err` in the work dir
  for a `gsnap FAILED (rc=...)` line, and see Known Issues.

A `gsnap.stats.txt` of only **9 bytes**, or a missing one, means gsnap produced an
empty alignment — treat the reads as untrustworthy. A healthy one is ~255 bytes and
contains a table of SAM FLAG counts.

---

# Part II — Metatranscriptomics (Nextflow phase)

Part II takes Part I's non-host reads and, in one Nextflow run: co-assembles contigs
(SPAdes), runs a host-filter BLASTn, a full `core_nt` BLASTn and a DIAMOND `blastx`
against NR, assigns taxonomy (`PROCESS_BLAST`/`PROCESS_DIAMOND`, using a `taxonomy.sif`
container) and draws an alluvial plot. It uses the **`rnaq-part2`** env with
`-profile slurm,singularity`.

> **Scope:** the Nextflow phase ends at the alluvial plot. Combining nt+nr, curation,
> and Salmon quantification are separate SLURM/R scripts (see the pipeline's
> `Metatranscriptome-Technical-Notes.md`) — out of scope here.

> **What a toy-DB run proves:** the plumbing (every step runs, containers work, formats
> are right) and a **positive control** — the toy DBs contain a known virus (plus host
> sequences), so a correct run recovers the seeded virus and calls everything else "no
> hit." In the example that leaves ~3 virus contigs at the end. It is **not** real
> biology — any result for anything other than the seeded virus is an artifact of the toy
> databases. Point the database paths at the toy DBs or the full DBs; the steps are the
> same.

---

## Part II · Step 1 — Confirm Part I actually produced reads

**Do not just check that the folders exist.** A failed Part I leaves the
`unmapped_reads/Paired/<accession>/` directories in place but *empty*, and the first
Part II step then fails with `ValueError: No non-host reads found for <BioProject>` —
a long way from the real cause.

```bash
PUB=$(sed -n 's/^publish-dir:[[:space:]]*//p' "$RUN"/params.zf-example.yaml)
find "$PUB/unmapped_reads" -name '*.fastq.gz' -size +132c -printf '%10s  %p\n' | sort -k2
```

Expect two non-empty `.fastq.gz` per sample (four total for the two-sample example).
If that prints nothing, or prints files of only a few bytes, **stop and fix Part I** —
see Part I Step 8 and Known Issues → "Part I reports success but produces no non-host
reads".

---

## Part II · Step 2 — Build the `unmapped_links` directory

Part II reads non-host reads via `--unmerged_accessions`, which the script points at
`$PROJECT_DIR/unmapped_links/` (your **run area** — the same `$RUN` you've used
throughout, *not* the repo). Build that directory and fill it with **symlinks to the
per-sample folders** from Part I's `unmapped_reads/Paired/`. Replace the target below
with your real Part I `publish-dir`:

```bash
mkdir -p "$RUN"/unmapped_links
cd "$RUN"/unmapped_links
# one `ln -s` per sample, pointing at the per-sample FOLDER under Part I's Paired output.
# The trailing "." links it into the current dir (unmapped_links/) keeping the accession name:
ln -s /path/to/your/RNAquarium_outputs/nonhost_zf-example/unmapped_reads/Paired/ERR2723539 .
ln -s /path/to/your/RNAquarium_outputs/nonhost_zf-example/unmapped_reads/Paired/ERR2723409 .
```

**Verify the links resolve** — a broken symlink looks fine to `ln` but fails at run
time, and is the #1 cause of an immediate Part II failure:

```bash
ls -laL "$RUN"/unmapped_links/    # -L follows links; each must show the target's real
                                  # contents, NOT "No such file or directory"
```

Each entry should read `ERR2723539 -> /…/unmapped_reads/Paired/ERR2723539/` and resolve.
If you see "No such file or directory," the `ln` target is wrong — check it against your
Part I `params.zf-example.yaml` `publish-dir:` value and redo.

---

## Part II · Step 3 — The assembly-group mapping

`zf-example-bioproject-mapping.json` maps each assembly group to its runs (samples in a
group are co-assembled). It's ready as-is for the example. Group IDs must match `PRJ` +
two uppercase letters + digits (the code validates this).

---

## Part II · Step 4 — Reference databases and the test config  (toy DBs: minutes; `apptainer build`: ~10 min)

**Toy route (fast plumbing + positive control)** — read `part2-test-databases.md`,
then build (uses the `rnaq-dbbuild` env from Setup):

```bash
conda activate rnaq-dbbuild
# edit the top of your copied build-toy-databases.sh: DB_DIR (+ optionally the virus
# accessions), then run on a LOGIN NODE (needs internet):
bash "$RUN"/build-toy-databases.sh
```

> If you see `mkdir: cannot create directory '/path': Permission denied`, you didn't edit
> `DB_DIR` — it's still the placeholder `/path/to/...`. Set it to a real scratch path.

It builds, in one `DB_DIR`: `host_toy.*`, `core_nt.*`, `nr.dmnd`, `nameNode.sqlite`, and
`taxdb.btd`/`taxdb.bti`. All databases go in one folder because BLAST finds its taxonomy
files by `cd`-ing into the mounted DB dir. **Only the seeded organisms are reportable.**

> The builder downloads **current** NCBI dumps for everything **except**
> `nameNode.sqlite`, which it builds from a pinned **Feb-2025** taxonomy dump. A current
> dump makes the public taxonomy script miscategorize viral-realm lineages as
> `other_Eukaryota` instead of `viruses`; Feb-2025 (the 75k run's vintage) fixes it and
> makes the run reproducible.

**Full route (real analysis)** — build the full databases per
[Inputs & Databases](https://czbiohub-sf.github.io/RNAquarium/inputs-and-databases.html)
and put `taxdb.btd`/`taxdb.bti` in each BLAST db directory.

Database parameters (all **snake_case**, matching the public `run_75k.sh`):

| Parameter | Points at |
| --------- | --------- |
| `--nt_dir` + `--nt_db_name` | host-filter db (`host_toy` for toy route) |
| `--ntfull_dir` + `--ntfull_db_name` | full `core_nt` |
| `--nr_dir` | DIAMOND NR (basename `nr` → `nr.dmnd`) |
| `--taxonomy_db` | taxonomizr `nameNode.sqlite` |

> Don't pass `--tax_btd`/`--tax_bti`/`--tax4blast`/`--taxonomy_db_fallback` — no process
> reads them; `taxdb.*` just needs to sit in the db directory.

**Config setup (both routes) — do once in your repo:**

```bash
conda deactivate            # build the container from a clean shell (not the dbbuild env)
cd "$REPO"/src/metatranscriptome
apptainer build containers/taxonomy.sif containers/taxonomy.def   # taxonomy steps need this
cp "$REPO"/docs/walkthrough/modules.walkthrough.config conf/modules.config  # test-scale cpu/mem/time
# revert for a real run: git checkout -- conf/modules.config
```

---

## Part II · Step 5 — Edit `run-part2-metatranscriptome.sh`

This script uses `#SBATCH --chdir` (it must start in the pipeline dir so `main.nf` and
`./containers/taxonomy.sif` resolve), so edit **three literal `#SBATCH` paths** at the
top (SLURM headers can't use variables):

- `--chdir=/path/to/RNAquarium/src/metatranscriptome`
- `-o /path/to/RNAquarium/metatranscriptome_slurm.out/slurm-%j.out` (absolute)
- `-e /path/to/RNAquarium/metatranscriptome_slurm.out/slurm-%j.err` (absolute)
- `--partition=` — your cluster's partition.

Then the body variables: `PROJECT_DIR`, `UNMERGED_ACC` (the `unmapped_links/` dir),
`BIOPROJ_MAP` (the JSON), and `DB_DIR` + the `NT_*`/`NTFULL_*`/`NR_*`/`TAXONOMY_DB`
paths (for the toy route, `NT_DB_NAME` is `host_toy`).

**SITE SETUP block:** set `module load anaconda` to your conda module (for `rnaq-part2`).
For the **container runtime**, check `command -v apptainer` / `command -v singularity` —
if either prints a path it's on PATH (leave the `module load` line commented); if
neither, uncomment and set it to your cluster's module name. The script uses
`-profile slurm,singularity`; if your cluster has apptainer but no `singularity` command,
change it to `slurm,apptainer` (see Known Issues).

---

## Part II · Step 6 — Run via `sbatch`  (~5–15 min)

The log dir must exist before submit (absolute `-o`/`-e` paths):

```bash
conda activate rnaq-part2                                   # provides Nextflow 24.10.5
mkdir -p /path/to/RNAquarium/metatranscriptome_slurm.out   # edit to your repo path
sbatch "$RUN"/run-part2-metatranscriptome.sh
squeue --me
```

`--chdir` handles the working directory, so you can submit from anywhere. (The
`conda activate rnaq-part2` matches the other run steps for consistency; the script's
SITE SETUP block also activates it inside the job, so the job is covered either way.)

---

## Part II · Step 7 — Find your outputs (and collect the unpublished ones)

`publish_dir` gets the SPAdes contigs (`single_end/`, `paired_end/`), non-host contig
FASTAs (`non_zf_hum_fa/`), the search tables (`nt_blast/`, `nr_diamond/`), and the HTML
run reports.

### The taxonomy tables and alluvial plot are NOT published

On the public version `PROCESS_BLAST`, `PROCESS_DIAMOND` and `ALLUVIAL_PLOT` have no
`publishDir`, so their outputs — including the headline result files — stay in the
Nextflow **work dir** and never reach `publish_dir`. This is expected, not a failure.

**Collect them into your run area with one block.** It reads the task work dirs out of
the Nextflow log, so you never have to hunt for hash directories by hand:

```bash
RUN=/path/to/your/run-area                     # same $RUN you've used throughout
LOG=$(ls -t "$RUN"/logs/nextflow_part2_zf-example.log* | head -1)
DEST="$RUN/part2_taxonomy_outputs"
echo "log:  $LOG"

# 1. Which unpublished tasks completed, and where did they run?
sed -n 's/.*name: \([^;]*\); status: COMPLETED;.*workDir: \([^ ]*\) started.*/\1\t\2/p' "$LOG" \
  | grep -E 'PROCESS_BLAST|PROCESS_DIAMOND|ALLUVIAL_PLOT|TREEMAP'

# 2. Copy their result files into $DEST, one subfolder per task
mkdir -p "$DEST"
while IFS=$'\t' read -r name wd; do
    sub="$DEST/$(echo "$name" | sed 's/[^A-Za-z0-9._-]/_/g')"
    mkdir -p "$sub"
    find "$wd" -maxdepth 1 -type f \
         ! -name '.command.*' ! -name '.exitcode' ! -name 'versions.yml' \
         -exec cp -p {} "$sub"/ \;
    echo "  $(ls -1 "$sub" | wc -l) file(s) -> $sub"
done < <(sed -n 's/.*name: \([^;]*\); status: COMPLETED;.*workDir: \([^ ]*\) started.*/\1\t\2/p' "$LOG" \
         | grep -E 'PROCESS_BLAST|PROCESS_DIAMOND|ALLUVIAL_PLOT|TREEMAP')

# 3. What you ended up with
find "$DEST" -type f -printf '%10s  %p\n' | sort -k2
```

The `status: COMPLETED` filter matters — it skips failed and retried attempts, so you
only collect work dirs that actually produced output.

**If the log is gone**, the same files can be found by name instead:

```bash
find "$RUN/scratch/rnaq-part2-work-zf-example" \
     \( -name 'taxonomy_hits_*' -o -name '*alluvial*' -o -name '*treemap*' \) \
     -printf '%10s  %p\n' | sort -k2
```

### What you should see

```
taxonomy_hits_nonhost_alluvialplot_all.png / .pdf   <- alluvial plot
taxonomy_hits_nonhost_treemap.png / .pdf            <- treemap
taxonomy_hits_viruses*.tsv                          <- the seeded virus
taxonomy_hits_chordates.tsv                         <- host / Danio hits
taxonomy_hits_<category>.tsv                        <- one per taxonomic category
```

On a toy-database run, **`taxonomy_hits_viruses*.tsv` should be non-empty** — that is
the positive control passing, meaning the seeded virus was recovered. A quick look:

```bash
head -3 "$DEST"/*ALLUVIAL_PLOT*/taxonomy_hits_viruses*.tsv
```

Every `taxonomy_hits_*.tsv` coming back empty is a different problem — see Known Issues
("Part II: taxonomy tables come back empty").

---

---

# Troubleshooting & Known Issues

Most Part I friction is handled automatically by the walkthrough's config/stats swaps
(Part I Step 6) — the `--qos`/`--tmp` strip, conda channels, the gsnap cpus guard, the
`.fna`→`.fa.gz` rename, the `fastq-numfilter` install, and the end-of-run stats-summary
fix are all baked in. The entries below are the ones a **correct run on a different
cluster can still legitimately hit** (mostly site-specific).

## Quick reference

| Symptom | Fix |
| ------- | --- |
| `--mkdir` error building Part I envs | ensure `mamba=1.5.*` in `rnaq-part1` (not `0.1.2`); run `-profile slurm,mamba` |
| `conda: command not found` inside the job | set your Case A/B/C line in the script's SITE SETUP block |
| `install-deps.nf` fails: env exists | remove the leftover env — if `conda env list` shows it **unnamed**, `-n` cannot work; use `conda env remove -p <path>`. See below |
| Part II: `Lmod ... unknown module "singularity"` | it's likely already on PATH — comment the module line; or use your cluster's name / `-profile apptainer` |
| **Part II: `No non-host reads found for <BioProject>`** | **Part I produced nothing; the cause is in Part I, not Part II — see below** |
| Part I says `Goodbye` but `unmapped_reads/` is empty | gsnap failed and was hidden by `errorStrategy 'ignore'` — see below |
| `gsnap: command not found` when testing by hand | Part I tools live in `src/nonhost/bin/`, not the conda env — see below |
| `gsnap.nosimd` fails in the self-test | real upstream GMAP bug, but only reachable on pre-2008 CPUs — run `bin/cpuid` to confirm yours has AVX2/SSE4.2 |
| Part II: taxonomy tables empty | search DB lacks taxonomy metadata — see below |
| Part II: taxonomy/alluvial missing from `publish_dir` | expected on public version — read them from the work dir (see Part II Step 7) |

## Known issues, explained

### Conda not available inside the batch job

A batch job starts a fresh shell that doesn't inherit your interactive conda. In the
script's SITE SETUP block, set the **make-conda-available** line for your case
(`module load <name>` A / `source ~/…/conda.sh` B / nothing C), then the always-run
`source "$(conda info --base)/…/conda.sh"` + `conda activate` lines. On many clusters A
plus those lines are both used — that's normal, not a contradiction.

### `mamba=1.5.*` (not `0.1.2`)

Part I's Nextflow 23.10.1 calls `conda create --mkdir`, which modern conda (≥24)
dropped; the older `mamba 1.5.x` resolver handles it. Pinning `mamba=1.5.*` (done in the
Step-setup env) avoids both the `--mkdir` error and a bogus `mamba 0.1.2`.

### `install-deps.nf` fails: `prefix already exists: .../seq-detective`

A leftover `seq-detective` conda env from a previous attempt blocks the rebuild. Remove
it and re-run — but **check `conda env list` first**, because the obvious command often
cannot work:

```bash
conda env list
```

```
                       /path/to/anaconda/23.1.0-3/x86_64/envs/seq-detective   <- NO NAME
rnaq-part1             /path/to/anaconda/25.3.1/x86_64/envs/rnaq-part1
```

**A blank name means the env is outside your active anaconda root's `envs_dirs`.** Conda
only assigns names to envs inside the root it is currently configured for, so for a
blank-name entry `conda env remove -n seq-detective` has no name to match and will
always fail — it is not a matter of it "sometimes not finding it". Note the root
versions differ above (`23.1.0-3` vs `25.3.1`): that is the tell. It happens because
`module load` (of anaconda, or of mamba during the build step) can switch roots between
sessions, so the env was created under one root and you are now using another.

Remove it by **path**, using the path from `conda env list` or from the error message:

```bash
conda env remove -p /path/to/anaconda/23.1.0-3/x86_64/envs/seq-detective
conda env list | grep -c seq-detective        # expect 0
nextflow run install-deps.nf
```

If conda still refuses, the env directory can be deleted directly — but **look before
you remove**, and be certain the path is an env and not something else:

```bash
ls -la /path/to/anaconda/23.1.0-3/x86_64/envs/seq-detective   # conda-meta/ should be here
```

The same reasoning applies to any unnamed entry in `conda env list`: Nextflow's own
per-process envs (`.../conda/env-<hash>`) are unnamed for exactly this reason, and they
belong to a run area rather than to you — leave those alone unless you are cleaning up
that specific run.

### `fastq-numfilter` missing from `bin/`

`setup-minimal.sh` builds it but (in some versions) doesn't copy it. Copy it in:
`find work -name fastq-numfilter -type f -exec cp {} bin/ \; -quit` (or use the newer
`setup-minimal.sh`).

### Part II: `Lmod ... unknown module "singularity"`

Your cluster has no module by that name — but Apptainer is usually already on PATH.
Check `command -v apptainer` / `command -v singularity`; if either prints a path, keep
the `module load` line commented. If only `apptainer` exists (no `singularity`
command), also change the profile to `-profile slurm,apptainer`.

### Part II: taxonomy tables come back empty (toy-DB runs)

Empty taxonomy output (no rows, no error) means the search DB lacks taxonomy metadata:
the scientific-name column comes back `N/A` and the R post-processing filters every row.
Ensure `taxdb.btd`/`taxdb.bti` sit **in the same directory** as the BLAST db
(`--ntfull_dir`/`--nt_dir`), and that `nr.dmnd` was built with
`diamond makedb --taxonmap … --taxonnodes … --taxonnames …`. `build-toy-databases.sh`
handles both.

### Part II: taxonomy tables / alluvial plot missing from `publish_dir`

Expected on the public version — the taxonomy/plot processes have no `publishDir`, so
their outputs stay in the Nextflow work dir. See Part II Step 7 for how to find them.
(A toy run's `taxonomy_hits_viruses*.tsv` being non-empty there is the
positive control passing.)

### Building the toy DBs: conda R / `taxonomizr` won't load

`build-toy-databases.sh` builds `nameNode.sqlite` with an `Rscript` that does
`library(taxonomizr)`, from the `rnaq-dbbuild` conda env (`r-taxonomizr`). If conda's R
misbehaves on your cluster (a known source of friction on some sites), use your own R
instead: `module load <your R module>` (e.g. `r/4.4`), then once
`R -e 'install.packages("taxonomizr", repos="https://cloud.r-project.org")'`. The
builder only needs an `Rscript` on `PATH` that can `library(taxonomizr)` — from conda or
a module.

### Part I reports success but produces no non-host reads

**Symptom.** Part I prints `Execution complete -- Goodbye` and exits 0. The
`unmapped_reads/Paired/<accession>/` directories exist but are empty or nearly so, and
`gsnap.stats.txt` is 9 bytes instead of ~255. Part II then fails at its first step with
`ValueError: No non-host reads found for <BioProject>`.

**Why it's silent.** Two things compound:

1. `withName: gsnap` resolves `errorStrategy` to `'ignore'` for any non-exit-9 failure,
   and `main.nf` deliberately routes around a failed gsnap. So a failed gsnap does not
   fail the workflow.
2. `step.7.gsnap.nf` guards its own fallback with `if [[ $? > 0 ]]` placed *after* a
   `set -e`. `$?` is therefore `set -e`'s status, not gsnap's, and `[[ ]]` makes `>` a
   string comparison rather than numeric. The test is never true, so the
   `gsnapSkipped` path — which exists precisely to pass reads through unchanged when
   gsnap fails — can never execute.

**The fix** is the `step.7.gsnap.walkthrough.nf` swap in Part I Step 6b. It captures
`gsnap_rc=$?` immediately, compares with `-gt`, and additionally treats an empty or
missing SAM as failure. With it in place a gsnap failure produces
`...gsnapSkipped.fastq.gz` and a `gsnap FAILED (rc=...)` line in the step's
`.command.err`, instead of silence.

**To diagnose an occurrence**, read the trace from that run — it records the status and
exit code of every task:

```bash
TRACE=$(ls -t "$RUN"/reports/part1-*/trace-*.txt | head -1)
grep -i gsnap "$TRACE"
```

`status=COMPLETED, exit=0` means gsnap's *script* succeeded and the problem is in what
it wrote; `FAILED`/`ABORTED` means the task itself failed and was ignored. The `hash`
column gives the work directory to inspect.

### `gsnap.nosimd` is broken upstream (but almost certainly not your problem)

`gsnap` is a dispatcher that picks a CPU-specific binary at runtime. In GMAP 2025-04-19
the **`gsnap.nosimd`** variant is defective: on this walkthrough's own example data it
exits 9 and writes no SAM, dying on a single short read —

```
Problem sequence: ERR2723539.308 (40 bp)
```

— while `gsnap.avx2` and `gsnap.sse42` align the identical input correctly (6106
records, 5788 unmapped). Reproduced on two independently compiled builds, so it is an
upstream defect, not a bad local compile.

**How likely are you to hit it? Very unlikely.** The dispatcher only falls back to
`nosimd` when the CPU has **neither AVX2 nor SSE4.2** — that means hardware older than
2008. Any x86-64 cluster node in service today supports at least SSE4.2, and anything
from 2013 onward has AVX2. This is recorded for completeness, not as an expected
failure.

**Check in five seconds** — GMAP ships its own detector, and `setup-minimal.sh` installs
it alongside the aligners:

```bash
"$REPO"/src/nonhost/bin/cpuid                          # login node
srun -p <your-partition> "$REPO"/src/nonhost/bin/cpuid # a compute node
```

If that reports AVX2 (or SSE4.2), the broken variant will never be selected and you can
ignore this entry. You can also confirm after the fact: the gsnap step's `.command.log`
records which binary actually ran, e.g. `called with args: gsnap.avx2 ...`.

**If you genuinely are on such hardware**, gsnap fails and — without the
`step.7.gsnap.nf` swap from Part I Step 6b — does so silently. With the swap you get
`...gsnapSkipped.fastq.gz` files and a `gsnap FAILED (rc=9)` line, and the run stays
usable.

### `gsnap: command not found` (and the other bundled tools)

Part I's aligners and helpers are **not in the `rnaq-part1` conda env**. The config
deliberately avoids the gmap containers, so `setup-minimal.sh` compiles them into
`<repo>/src/nonhost/bin/`, and they are only on `PATH` because the submission script
does:

```bash
export PATH="$PATH:${REPODIR}/src/nonhost/bin"
```

Activating `rnaq-part1` alone will not give you `gsnap`. To use them interactively:

```bash
export PATH="$PATH:$REPO/src/nonhost/bin"
command -v gsnap fastq-namefilter seq-detective
```

### First run needs outbound network for the `nf-validation` plugin

Both parts declare `plugins { id 'nf-validation' }`. On the **first** run Nextflow
downloads it into `~/.nextflow/plugins/` — so a cluster whose login and compute nodes
have no outbound internet will fail at startup, before any pipeline step runs.

Once cached it is never fetched again, which means this is easy to miss when testing on
a machine that has run Nextflow before. To check whether you are relying on a cache:

```bash
ls ~/.nextflow/plugins/          # nf-validation-*  = already cached
```

If your compute nodes are offline but a login node is not, running Nextflow once from
the login node populates the cache for later batch jobs. Otherwise the plugin directory
can be copied in from a machine that does have access.

### Nextflow options differ between versions

Part I uses Nextflow **23.10.1** and Part II uses **24.10.5**, and their CLI option
sets are not the same — nor the same as the current online documentation. For example
`nextflow config` accepts `-params-file` on newer releases but reports
`Unknown option: -params-file` on 23.10.1. Treat `nextflow help <command>` on the
installed build as authoritative. Note also that `-c` and `-log` are **launcher**
options and belong before the subcommand (`nextflow -c my.config run ...`).

### `--nice` rejected (rare)

The Part I test config keeps `--nice=25` on the `download` step. If your cluster rejects
it, edit that `withName: download` line in your copied `nextflow.config` to
`clusterOptions = ''`.

---

# Roughly how long does it take?

Ballparks on a SLURM cluster with the tiny example (varies with load/network):

| Step | Rough time |
| ---- | ---------- |
| Part I — first-run aligner index build (one-time; cached after) | **~25 min** |
| Part I — per-sample steps after indexes | a few min |
| Part II — `apptainer build taxonomy.sif` (one-time) | **~10 min** |
| Part II — full Nextflow run (toy DBs) | **~5–15 min** |

Everything else (env creation, `setup-minimal.sh` ~10 min, `install-deps.nf` ~5 min,
downloads) is a one-time setup cost. The toy-DB build itself is minutes.

**Budget 1 to 1.5 hours end to end for a first pass.** The table above is machine time;
the rest is you. Every template you copy into your run area has paths, a partition name
and a genome size to fill in, and each one rewards being read before it is run —
`params.zf-example.yaml`, both submission scripts, and (for a toy-database run)
`build-toy-databases.sh`. Reading them is not busywork: most of the failure modes in
Known Issues announce themselves as a placeholder you didn't replace or a site setting
that doesn't match your cluster, and catching one at edit time costs a minute where
catching it mid-run costs a rebuild. A second pass, with the paths already known and the
indexes cached, is far quicker.
