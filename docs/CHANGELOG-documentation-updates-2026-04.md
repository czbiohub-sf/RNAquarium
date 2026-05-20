# Documentation Updates — April 2026

Summary of all documentation changes made during the Sphinx documentation overhaul and audit-driven update session.

---

## New files created

### `docs/inputs-and-databases.md`
New consolidated page documenting all input files and external database requirements for both pipeline parts. Covers host genome (with GRCz12tu as current reference), GTF annotations, ERCC spike-ins, SRA accession lists, contaminant sequences, pre-built aligner indexes, NCBI NT/NR databases, Taxonomizr setup, BLAST taxonomy files, and BBDuk adapters. Includes download URLs, format requirements, version recommendations, a summary table with approximate sizes, and a disk space planning section. Added to the Sphinx sidebar as the first entry under Part I.

### `docs/nonhost-readme.md`
New Sphinx-integrated copy of `src/nonhost/README.md`, added to the Part I toctree. Subsequently updated with all the changes described below.

### `docs/metatranscriptome-readme.md`
New Sphinx-integrated copy of `src/metatranscriptome/README.md`, added to the Part II toctree.

### `docs/AUDIT-detailed-75k-comparison.md`
Detailed audit report comparing the actual 75k production run code (`actual75kcode/`) against the GitHub repo (`src/`). Documents genome reference update path (GRCz11 → GRCz12tu via params override), tool version mismatches, undocumented features, and prioritized update recommendations.

### `docs/_images/RNAquarium-Pipeline-for-Repo_v042026.png`
New pipeline diagram image replacing the previous SVG, showing both Workflow I (Transcriptomic + Filtering) and Workflow II (Metatranscriptomic) with updated step names and tool versions.

---

## Files rewritten

### `docs/pipeline-explanation.md`
Complete rewrite. Removed the `[OUTDATED]` tag. Now accurately documents the current pipeline:

- **Step 0:** Index generation for 5 aligners (HISAT2, STAR, STAR+ERCC, Bowtie2, GSNAP) plus optional kallisto contamination index
- **Step 1:** SRA download with seq-detective barcode detection and filtering (previously undocumented)
- **Step 2:** fastp adapter trimming and quality filtering (PRICE removed — was documented but no longer in code)
- **Step 2.5:** Kallisto contamination filter via pseudoalignment (previously undocumented)
- **Steps 3–7:** HISAT2 → STAR → Bowtie2 → czid-dedup → GSNAP sequential host filtering
- **Step 3b:** Parallel host gene counting branch with featureCounts (default) or HTSeq, plus optional CRAM output (previously undocumented)
- Added output directory structure documentation

### `docs/parameters.md`
Complete rewrite. Reorganized into logical sections (Reference genome, Pre-built indexes, Input/output, Pipeline flow control, Resource/performance, Filtering, Publishing, Advanced). Removed ~10 stale parameters that no longer exist in code (`priceOptions`, `fastpOptions`, `starOptions`, `starIndexGenOptions`, `hisatIndexGenOptions`, `hisatOptions`, `sraPrefetchOptions`, `fastqDumpOptions`, `samtoolsSortOptions`, `starCountOptions`, `publishPricefiltered`). Added ~20 active parameters that were undocumented (`cleanupIntermediate`, `contamFa`, `dedupMinLen`, `dedupPercentLen`, `extraAdapters`, `htseqCount`, `kbContamIndexes`, `kbRetainMixed`, `nxfUnstageHack`, `outputCram`, `publishBowtie`, `publishDedup`, `publishGsnap`, `publishKallisto`, `publishQCfiltered`, `retainMixed`, `sdFilterMates`, `seed`, `starSjdbOverhang`, `starThreadsLarge`, `starThreadsSmall`, `starUseSharedMem`, `backupScratchHack`, `maxCpus`, `maxMemory`). Added descriptions for all parameters. Added note explaining that defaults still reference GRCz11 while the 75k run used GRCz12tu via params file override.

### `docs/nonhost-readme.md` and `src/nonhost/README.md`
Rewrote dependency list: removed PRICE (no longer used in pipeline). Added 7 tools not previously documented (kallisto, bustools, kb-python, seqtk, subread/featureCounts, updated czid-dedup version). Dependency table now shows both conda and container version numbers. Updated all genome examples from GRCz11/Ensembl release 108 to GRCz12tu with genome-size `1448808683`. Updated pipeline invocation examples to use `-profile slurm,mamba` and `-params-file` pattern matching the actual 75k run. Updated feature list to mention kallisto filter, seq-detective, and featureCounts.

---

## Files with targeted edits

### `docs/index.md`
- Replaced pipeline diagram from `RNAquarium-Pipeline-for-Repo_v3.svg` to `RNAquarium-Pipeline-for-Repo_v042026.png`
- Split single toctree into two sections: "Part I: Transcriptomic + Filtering" and "Part II: Metatranscriptomics"
- Added `inputs-and-databases`, `nonhost-readme`, and `metatranscriptome-readme` to the toctrees
- Added the two metatranscriptome technical notes docs (previously missing from navigation)

### `README.md` (root)
- Updated pipeline diagram reference from old SVG to new PNG

### `docs/conf.py`
- Updated project name to "RNAquarium"
- Updated copyright to "2026, Biohub"
- Updated author list to all 11 contributors
- Added RTD theme options (navigation depth, sticky nav, non-collapsing sidebar)
- Added `colon_fence` MyST extension and heading anchors
- Configured logo as `_static/logo.png`

### `docs/gettingstarted.md`
- Updated `genome-size` example from `1396431182` (GRCz11) to `1448808683` (GRCz12tu)
- Updated example Nextflow process list: removed `priceseqfilter`, `prefetch`, `fastq_dump`; added `kb_negative`, `feature_count`
- Fixed two broken wiki-style links (`[[params.html]]` → proper Markdown links)

### `docs/tips.md`
- Updated genome filename example to show both GRCz11 and GRCz12tu naming patterns

### `docs/metatranscriptome-readme.md`
- Updated Part I reference from "host mapping pipeline" to "transcriptomic + filtering pipeline"

---

## Naming convention change

"Host Mapping" → **"Transcriptomic + Filtering"** for Part I throughout all documentation. Updated in: index.md sidebar caption, pipeline-explanation.md title, parameters.md title, inputs-and-databases.md heading and intro, metatranscriptome-readme.md features list.

---

## Not changed (deliberate)

- **`docs/technical-notes.md`** — Content still current (Nextflow unstaging workarounds, dependency notes). No PRICE references found.
- **`docs/Metatranscriptome-Technical-Notes.md`** — Content verified as accurate, including the Taxonomizr feb2025 → august2025 timeline.
- **`docs/Metatranscriptome_virus-steps-Technical-Notes.md`** — Content verified as accurate.
- **`src/nonhost/setup.sh`** — Still downloads old tool versions (sra-tools 3.0.0, fastp 0.23.2, PRICE). This is a code change, deferred for separate discussion.
- **`docs/pipeline.pikchr`** — Orphan file, not referenced anywhere. Can be safely removed.
- **`docs/_images/RNAquarium-Pipeline-for-Repo_v3.svg`** and **`Q2_update_TLG3_v2_RNAquarium_pipeline_v1.svg`** — Old diagram files retained for historical reference.
