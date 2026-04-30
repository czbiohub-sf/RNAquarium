# RNAquarium Detailed Audit: actual75kcode vs GitHub Repo

**Date:** April 10, 2026
**Scope:** File-by-file comparison of `/actual75kcode/` against `/src/`, plus documentation update recommendations

---

## Executive Summary

The core pipeline code (main.nf, nextflow.config, all modules) is **identical** between the actual 75k run and the GitHub repo. The 75k run used `params.75k.3.yaml` to override defaults at runtime — this is how GRCz12tu replaced GRCz11 without changing the base code. The documentation drift is therefore not about code changes, but about:

1. **Default values in code/docs still reference GRCz11** — the 75k run overrode these, but the repo defaults are stale
2. **setup.sh downloads old tool versions** — the 75k run used conda/containers with newer versions
3. **New pipeline features (seq-detective, kallisto, featureCounts) are undocumented**
4. **Removed features (PRICE) are still documented**
5. **~20 new parameters are undocumented; ~10 documented parameters are gone**

---

## 1. The GRCz11 → GRCz12tu Story

### How the 75k run actually updated the genome

The code in `main.nf` still defaults to GRCz11:

```groovy
// main.nf lines 19-22 (STILL IN REPO)
params.refGenome = "Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz"
params.refGenomeGtf = "Danio_rerio.GRCz11.108.gtf"
```

The 75k run overrode this via the `-params-file` flag in `nextflow-submit-75k.sh`:

```bash
nextflow run main.nf -params-file params.75k.3.yaml -profile bruno,mamba ...
```

And `params.75k.3.yaml` specifies:

```yaml
ref-genome: .../GRCz12tu/GCF_049306965.1_GRCz12tu_genomic.fa.gz
ref-genome-gtf: .../GRCz12tu/GCF_049306965.1_GRCz12tu_genomic.gtf
genome-size: 1448808683  # was 1396431182 for GRCz11
```

### Evolution of params files across 75k attempts

| File | Key differences |
|------|----------------|
| `params.75k.yaml` (attempt 1) | GRCz12tu genome, `contam-fa: contaminant_danioRefseq_addgenePlasmid.fa`, `backup-scratch-hack: true`, parallel-downloads: 50, no pre-built indexes |
| `params.75k.2.yaml` (attempt 2) | Added pre-built index paths for all 5 mappers + kb-contam-indexes, `contam-fa: drerio_refseq_and_addgene_regions.fa` (updated), `backup-scratch-hack: false`, parallel-downloads: 120 |
| `params.75k.3.yaml` (attempt 3, final) | Same as 2 but parallel-downloads: 40, publish-dir changed to `host_mapping2/` |

### What needs to change in the repo

**14 files** still reference GRCz11. The recommended approach:

**Option A (Minimal):** Update `params.example.yaml` and docs to show GRCz12tu examples, but keep `main.nf` defaults generic (users must provide their own genome). This matches how the 75k run actually worked.

**Option B (Full update):** Change all defaults to GRCz12tu. This makes the repo zebrafish-specific by default.

Either way, these need updating:

| File | Change needed |
|------|--------------|
| `src/nonhost/params.example.yaml` | Update genome filename examples to GRCz12tu |
| `docs/nonhost-readme.md` | Update Ensembl example commands |
| `docs/parameters.md` | Update default values table |
| `docs/tips.md` | Update index directory naming examples |
| `src/nonhost/nextflow_schema.json` | Update schema defaults |

**New value:** `genome-size: 1448808683` (was `1396431182`)

---

## 2. Tool Version Mismatches: setup.sh vs Conda/Container Reality

The `setup.sh` script downloads tools manually for environments without conda/containers. It is identical in `actual75kcode` and the repo, but its versions are **far behind** what the pipeline actually uses via conda:

| Tool | setup.sh version | Container version | Conda version | Documented version |
|------|-----------------|-------------------|---------------|-------------------|
| **sra-tools** | 3.0.0 | 3.2.1 | 3.2.1 | 3.0.0 |
| **fastp** | 0.23.2 | 0.23.4 | 0.24.1 | 0.23.2 |
| **STAR** | 2.7.10b | 2.7.10b | 2.7.10b | 2.7.10 |
| **samtools** | 1.16.1 | 1.19.2 | 1.21 | 1.16.1 |
| **bowtie2** | *(not in setup.sh)* | 2.4.5 | 2.5.4 | 2.4.5 |
| **fastq-lengths** | 0.1.1 | — | — | *(not documented)* |
| **PRICE** | 140408 | — | — | Listed as dependency |

**Key observations:**

- `setup.sh` still downloads **PRICE**, which is no longer used in the pipeline
- `setup.sh` still downloads **sra-tools 3.0.0** but conda uses **3.2.1**
- fastp has **three different versions**: 0.23.2 (setup.sh), 0.23.4 (container), 0.24.1 (conda)
- samtools has **three different versions**: 1.16.1 (setup.sh), 1.19.2 (container), 1.21 (conda)
- The conda profile bundles many tools in a single environment spec (sra-tools, seqtk, jq, polars, fastp, hisat2, subread, samtools) — this is the de facto specification

**Recommendation:** Update `setup.sh` to match conda versions. Remove PRICE download. Add missing tools (seqtk, kallisto, bustools, kb-python, subread). Update `nonhost-readme.md` and `docs/parameters.md` to list the conda versions as canonical.

---

## 3. Undocumented Pipeline Features

### 3a. seq-detective (Step 1)

**In code:** `step.1.nf` — barcode detection process that identifies and filters technical reads and barcode sequences from scRNA-seq libraries.

**In docs:** Not mentioned anywhere.

**Recommendation:** Add to pipeline-explanation.md (or its replacement) and gettingstarted.md.

### 3b. Kallisto contamination filter (Step 2.5 / Step 8)

**In code:** `step.8.kallisto.nf` — uses kallisto pseudoalignment for host/contaminant filtering. Active in `main.nf` workflow. Requires `kb-contam-indexes` parameter.

Tools: kallisto 0.51.1, bustools 0.44.1, kb-python 0.29.1

**In docs:** Only visible as parameter defaults in `nextflow_schema.json`. No explanation of what it does or when to use it.

**Recommendation:** Add documentation explaining this as an additional host-filtering step using pseudoalignment.

### 3c. featureCounts / subread (Step 3b)

**In code:** `step.readcounts.nf` — subread 2.0.6 used for featureCounts as an alternative to HTSeq for host gene quantification.

**In docs:** Not mentioned. HTSeq is the only documented counting tool.

**Recommendation:** Document alongside HTSeq.

### 3d. Additional undocumented tools in conda spec

| Tool | Version | Purpose |
|------|---------|---------|
| seqtk | 1.5 | Read sampling/manipulation |
| polars | 1.16.* | Fast data processing |
| jq | (conda-forge) | JSON parsing |

---

## 4. Removed Feature Still Documented: PRICE

PRICE (PriceSeqFilter) is referenced in:

- `docs/pipeline-explanation.md` — described as part of Step 2
- `docs/nonhost-readme.md` — listed as a required dependency
- `docs/technical-notes.md` — mentioned in passing
- `src/nonhost/setup.sh` — still downloaded and compiled

But PRICE is **not used in any active Nextflow module**. The `install-deps.nf` script has PRICE installation commented out.

**Recommendation:** Remove from all documentation. Remove from `setup.sh`. Note in changelog that fastp replaced PRICE for quality filtering.

---

## 5. Parameter Documentation Gaps

### Parameters in code but NOT documented (~20)

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `cleanupIntermediate` | true | Remove intermediate files after completion |
| `contamFa` | '' | Contaminant FASTA for kallisto filtering |
| `dedupMinLen` | 50 | Minimum read length for deduplication |
| `dedupPercentLen` | 0.9 | Percent length threshold for dedup |
| `extraAdapters` | '' | Additional adapter sequences file |
| `htseqCount` | false | Use HTSeq instead of featureCounts |
| `kbContamIndexes` | '' | Pre-built kallisto contaminant index |
| `kbRetainMixed` | false | Keep mixed reads in kallisto filter |
| `nxfUnstageHack` | false | Workaround for Nextflow unstaging bug |
| `publishBowtie` | false | Publish Bowtie2 intermediate files |
| `publishDedup` | false | Publish dedup intermediate files |
| `publishGsnap` | false | Publish GSNAP intermediate files |
| `publishKallisto` | false | Publish kallisto intermediate files |
| `publishQCfiltered` | false | Publish QC-filtered files |
| `retainMixed` | false | Retain mixed single/paired reads |
| `sdFilterMates` | true | seq-detective mate filtering |
| `seed` | 0 | Random seed for reproducibility |
| `starSjdbOverhang` | 100 | STAR splice junction overhang |
| `starThreadsLarge` | 12 | STAR threads for large files |
| `starThreadsSmall` | 4 | STAR threads for small files |
| `starUseSharedMem` | false | STAR shared memory genome mode |

### Parameters documented but NOT in code (~10, stale)

| Parameter | Notes |
|-----------|-------|
| `fastpOptions` | Replaced by hardcoded fastp flags |
| `fastqDumpOptions` | Removed |
| `hisatIndexGenOptions` | Removed |
| `hisatOptions` | Removed |
| `priceOptions` | PRICE was removed entirely |
| `samtoolsSortOptions` | Removed |
| `starCountOptions` | Removed |
| `starIndexGenOptions` | Removed |
| `starOptions` | Removed |
| `sraPrefetchOptions` | Removed |

---

## 6. Part II (Metatranscriptome) Findings

### 6a. Core code is identical

All core files (main.nf, nextflow.config, conf/modules.config, all modules, all subworkflows, all bin/ scripts) are **byte-for-byte identical** between `actual75kcode/metatranscriptome` and `src/metatranscriptome`.

### 6b. Actual75kcode is a subset of repo

The repo has **71 files** vs **34 files** in actual75kcode. The missing files are:

- `containers/` directory (taxonomy.def Singularity definition)
- `scripts/post_processing/` directory (38 SLURM/R scripts for taxonomy, virus curation, quantification)

These post-processing scripts exist only in the repo because they were run separately on the HPC, not via Nextflow.

### 6c. Taxonomizr database version timeline (CONFIRMED)

The documentation at `Metatranscriptome-Technical-Notes.md` lines 429-435 **does correctly explain** the two-phase approach:

> "Initial runs used a February 2025 Taxonomizr database."
> "Fix: Download updated Taxonomizr database (August 2025) and re-run Steps 0a, 0b, 1, 2, etc."

This is accurate: the Nextflow pipeline used feb2025, and the standalone SLURM scripts were re-run with august2025. The previous audit incorrectly flagged this as a discrepancy — the documentation is correct.

### 6d. Minor differences in actual75kcode

| File | Difference |
|------|-----------|
| `scripts/create_json_from_mapping.py` | Paths adjusted for 75k dataset; some fields commented out |
| `modules/preblast1.nf` | Extra file in actual75kcode only (unused duplicate, not in repo) |
| `fastp_adapters_with9added.fasta` | Different directory location (root vs bin/) |

### 6e. Still-relevant Part II documentation gaps

- Diamond parameters (`blastx --ultra-sensitive --top 3`) not documented
- Host-only BLAST evalue (`1e-20`) and max_target_seqs (`1`) not documented
- R package versions in taxonomy container not pinned

---

## 7. Files in actual75kcode NOT in the Repo (Candidates for Addition)

### Part I nonhost

| File | Description | Add to repo? |
|------|-------------|-------------|
| `nextflow-submit-75k.sh` | SLURM submission script for 75k run | Yes — as a template/example |
| `params.75k.3.yaml` | Final 75k params with GRCz12tu | Yes — as `params.example.75k.yaml` |
| `params.75k.yaml` | First attempt params | Optional — historical reference |
| `params.75k.2.yaml` | Second attempt params | Optional — historical reference |
| `publish_download_stats.sh` | Extract download stats from NF logs | Optional — utility |
| `bin/` directory (65 files) | Pre-compiled binaries (gsnap, czid-dedup, etc.) | No — too large, use conda/containers |

### Part II metatranscriptome

| File | Description | Add to repo? |
|------|-------------|-------------|
| `modules/preblast1.nf` | Unused duplicate module | No |

---

## 8. Recommended Update Priority

### Tier 1: Critical for accuracy

1. **Update `params.example.yaml`** to show GRCz12tu filenames and genome-size `1448808683`
2. **Update documentation genome references** (nonhost-readme, parameters.md, tips.md) to GRCz12tu
3. **Remove PRICE** from docs and setup.sh
4. **Update dependency version table** in nonhost-readme.md to match conda specs
5. **Add `params.75k.3.yaml`** (or sanitized version) as an example params file

### Tier 2: Important for completeness

6. **Document seq-detective, kallisto filter, featureCounts** in pipeline docs
7. **Add ~20 missing parameters** to parameters.md
8. **Remove ~10 stale parameters** from parameters.md
9. **Update `setup.sh`** tool versions to match conda (or deprecate setup.sh in favor of conda)
10. **Rewrite pipeline-explanation.md** (currently marked OUTDATED)

### Tier 3: Nice to have

11. Document Diamond/BLAST parameter choices in metatranscriptome tech notes
12. Pin R package versions in taxonomy container
13. Add `nextflow-submit-75k.sh` as a submission template
14. Reconcile fastp version between container (0.23.4) and conda (0.24.1)
15. Update `nextflow_schema.json` defaults
