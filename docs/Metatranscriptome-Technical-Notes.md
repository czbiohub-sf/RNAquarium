# RNAquarium 75k Metatranscriptome – Technical Notes

These notes summarize the key technical decisions and workarounds needed to reproduce the 75k metatranscriptome run of the RNAquarium pipeline, focusing on Part II (metatranscriptome). References to internal HPC paths have been anonymized and may need adaptation on other systems.

---

## Overview

Part II of the RNAquarium pipeline ("metatranscriptome") consumes **unmapped reads** from the host-mapping stage (both **single-end** and **paired-end**) and produces:

- Non-host contigs assembled by **SPAdes**
- Host-filtered BLASTn results (Danio + human)
- Full **nt BLAST** and **Diamond nr** results
- Curated, BBDuk-filtered non-host contigs
- Virus-focused outputs (contigs, clusters, trees)
- Virus quantification via **Salmon**

**Scale:** ~75k SRA runs, organized into ~2000 Bioprojects.

**Strategy:**

- One main Nextflow run for most Bioprojects.
- A secondary "last9" Nextflow run for the 9 largest Bioprojects (further subdivided) plus 7 Bioprojects that were missed during a resume.
- Heavy use of Nextflow caching/resume and manual SLURM control for problematic steps.
- Post-search adapter/contaminant filtering using BBDuk at the **contig** level (not reads) to avoid a full re-run.

---

## Pipeline Architecture: Nextflow vs SLURM Scripts

**Important:** The Nextflow pipeline only covers the initial phases of the metatranscriptome workflow. All downstream processing is handled by standalone SLURM/R scripts located in `src/metatranscriptome/scripts/post-processing/`. These should be moved into your working output folder (the PUB_DIR in run_xxx.sh file), with a subset inside output subfolders.

### What Nextflow Handles

1. SPAdes assembly of non-host reads into contigs
2. Host-filtering BLASTn (Danio + human removal)
3. Full nt BLASTn search
4. Diamond nr search

### What SLURM Scripts Handle (Post-Nextflow)

| Step | SLURM Script | R Script(s) | Description |
|------|--------------|-------------|-------------|
| 0a | `slurm_0a_nt_blast_Rcommands.sh` | `nt_blast/nt_blast_processing.r` | Parse BLASTn outputs with Taxonomizr |
| 0b | `slurm_0b_nr_diamond_Rcommands.sh` | `nr_diamond/nr_diamond_processing.r` | Parse Diamond outputs with Taxonomizr |
| 1 | `slurm_1_combinecontigs.sh` | `1_combinecontigs_alluvial_treemap.r` | Combine NT + NR results, LCA calculations |
| 1b | `slurm_1b_combinecontigs.sh` | `1_combinecontigs_alluvial_treemap_bbdukfilter.r` | Combine with BBDuk filtering; also calls `taxname_to_taxid_*.sh` |
| 2 | `slurm_2_createfastas.sh` | `3_addsequence.r` | Create FASTAs for each broad group |
| 2 | `slurm_2_createfastas_separate.sh` | `3_addsequence_separatefastas.r` | Variant: separate FASTAs |
| 2b | `slurm_2b_bbdukflag.sh` | — | BBDuk adapter masking on contigs |
| 3 | `slurm_3_combinecontigsandsequence.sh` | `3_addsequence.r` | Add sequences to taxonomy tables |
| 4 | `slurm_4_counting_anddarkmatter.sh` | `4_alluvial_withdarkmatter.r` | Counting and dark matter analysis |
| 4b | `slurm_4b_counting_anddarkmatter.sh` | `4b_alluvial_withdarkmatter.r` | Counting and dark matter with BBDuk filter |
| 4t | `slurm_4t_general_metacoder_trees.sh` | `general_metacoder_trees.r` | General metacoder trees (fungi, etc.) |

Helper scripts called by step 1b:
- `taxname_to_taxid_nonhost_twoparts.sh` – adds `taxid_lca_NTorNR` for all non-host contigs
- `taxname_to_taxid_viruses0_twoparts.sh` – adds `taxid_lca_NTorNR` for virus contigs

Virus-specific scripts (see **Metatranscriptome_virus-steps-Technical-Notes.md**):

| Step | SLURM Script | R Script(s) | Description |
|------|--------------|-------------|-------------|
| v5A | `slurm_virus_5A_curation.sh` | `virus_curation_plots_addsequenceA.r` | Add realm taxonomy (if needed) |
| v5 | `slurm_virus_5_curation.sh` | `virus_curation_plots_addsequence.r` | Virus filtering and annotation |
| v6 | `slurm_virus_6_metacodertrees.sh` | `virus_metacoder_trees.r` | Virus metacoder trees |
| v7 | `slurm_virus_7_clustering.sh` | `virus_clustering_curation.r` | Virus clustering |
| v7b | `slurm_virus_7b_clustering_afterbbduktrim.sh` | `virus_clustering_curation.r` | Clustering with BBDuk filter |
| v8 | `slurm_virus_8_clusteringII.sh` | `virus_clustering_curation.r` | Clustering phase II |
| v9 | `slurm_virus_9_clustering_creatingfastas.sh` | `virus_clustering2_outputingmanyfastaclusters.r` | Create cluster FASTAs |
| v10 | `slurm_virus_10_clustering_minimap.sh` | — | Minimap clustering |

Scripts in `nr_diamond/` subfolder:
- `slurm_diamond_taxid_update_blastdbcmd.sh` – recover missing Diamond taxids
- `diamond_taxid_update_blastdbcmd.py` – Python script for taxid lookup
- `diamond_taxid_moveold_andrename.py` – rename updated files

Even with no errors, users should expect to run SLURM scripts sequentially after Nextflow completes.

---

## Inputs & Basic Setup

### Unmapped Read Inputs

All unmapped reads from Part I were symlinked into a single directory combining SE + PE runs:

```
<PIPELINE_ROOT>/src/metatranscriptome/all_unmapped_links
```

Original locations (for reference) are in the output folder of the Host Mapping (Part I) pipeline portion:

```
# Single-end
<OUTPUT_ROOT>/host_mapping/unmapped_reads/Single

# Paired-end
<OUTPUT_ROOT>/host_mapping/unmapped_reads/Paired
```

Both `run_75k.sh` and `create_json_from_mapping.py` must reference `all_unmapped_links`.

### SRA Metadata & Run→Bioproject Mapping

We rely on the NCBI SRA accessions master file and a custom mapping file:

```bash
# Download full SRA accessions table (once, global)
cd <SRA_METADATA_DIR>
wget ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab
```

Run-to-Bioproject mapping:

- Curate a run → Bioproject TSV: `RNaquarium_run_to_bioproject_mapping.tsv`
- Keep this mapping and associated scripts in `<SRA_METADATA_DIR>`

The mapping TSV is converted to JSON with:

```bash
module load data.science
module load anaconda
python <PIPELINE_ROOT>/src/metatranscriptome/scripts/create_json_from_mapping.py
```

The resulting JSON encodes the Bioproject/run relationships used by Nextflow.

### Containers

We use an apptainer container for taxonomy-related steps:

```bash
cd <PIPELINE_ROOT>/src/metatranscriptome
apptainer build containers/taxonomy.sif containers/taxonomy.def
```

### General Environment

- Nextflow (version pinned in `run_75k.sh`)
- Modules: `data.science`, `anaconda`, plus any local Nextflow/Java dependencies defined in your environment

---

## Nextflow Pipeline Specifics

### Main Run Script

Primary launch script for the main run:

```bash
cd <PIPELINE_ROOT>/src/metatranscriptome
sbatch run_75k.sh
```

**Important:** Both `run_75k.sh` and `create_json_from_mapping.py` must reference the same `all_unmapped_links` directory.

### Key Config and Module Files

- `nextflow.config`
- `modules/*.nf` (notably `assembly.nf`, `preblast.nf`, `blast.nf`)
- `subworkflows/align_nt.nf`
- `conf/modules.config`
- `run_75k.sh`

---

## Critical Configuration Changes (Required for 75k Run)

### Nextflow Version Override Fix

An early failure was caused by an `NXF_VER` override variable. Remove or correct any override that forces an incompatible Nextflow version. Ensure `run_75k.sh` uses a tested Nextflow version.

### SPAdes Retry / Exit-Code Handling

SPAdes will occasionally fail on very large Bioprojects. We:

- Changed the SPAdes "hard failure" exit code to 100 to allow Nextflow to treat it as ignorable after sufficient retries.
- Increased retry attempts and dynamic walltime for difficult SPAdes tasks.

Example in `modules/assembly.nf`:

```groovy
// Change SPADES error code in assembly module to 100
// so that we can classify it as ignorable after a threshold.
// (Previously exit code 9 was treated as hard failure)
exit_code = 100
```

Example in `conf/modules.config` for SPAdes:

```groovy
withName: ASSEMBLE_SPADES_* {
  maxRetries = 25
  time = {
    task.attempt == 1 ? 4.hour :
    task.attempt == 2 ? 4.hour :
    task.attempt == 3 ? 48.hour :
    80.hour
  }
  errorStrategy = {
    (task.exitStatus == 100 || task.attempt > 20) ? 'ignore' : 'retry'
  }
}
```

Key conventions:

- Exit code 100 is treated as "skip this task but continue the pipeline".
- We used this to skip the largest 9 bioprojects in the main run and handle them separately.
- To do so, we manually wrote `100` into `.exitcode` files in the relevant Nextflow work directories when necessary.

### Global Nextflow Retries

There was a global `maxRetries = 5` higher up in `modules.config` that caused premature failure. For 75k, we effectively allowed more SPAdes attempts via the block above and avoided early global failure by either removing or relaxing the global `maxRetries`.

### First BLAST Taxid Setting

The initial host BLAST step must restrict taxids to human + Danio. The bug we hit:

- `nextflow.config` had taxid parameters, but
- `modules/blast.nf` did not pass `-taxids`.

This must be fixed so that the first BLAST only returns host (Danio/human) hits for removal:

```bash
blastn \
  -db ${params.core_nt_db} \
  -query chunk.fasta \
  -taxids ${params.host_taxids} \
  ...
```

Missing `-taxids` caused spurious non-host hits in the host-filtering BLAST; the main run had to be stopped and resumed after adding the flag.

### Full BLAST Variable Wiring

Ensure `NTFULL_DIR` and `NTFULL_DB_NAME` are properly passed through:

- `subworkflows/align_nt.nf`
- `nextflow.config`
- `run_75k.sh` (e.g. `--ntfull_dir $NTFULL_DIR \` and `--ntfull_db_name $NTFULL_DB_NAME \`)

### Host-Only BLAST vs Full BLAST

- **First BLAST step:** Host-only (`-taxids` = human, Danio) to create non-host contigs.
- **Second BLAST step:** Full NT BLAST on non-hosts, using `core_nt` database.

---

## Handling the 9 Largest Bioprojects ("last9" Run)

### Motivation

With increased SPAdes retries and longer walltime, a small number of very large Bioprojects dominated runtime:

- 9 Bioprojects (~0.5% of Bioprojects) contained ~22% of total reads.
- All 9 were single-end and too big for SPAdes under default chunking.

### Approach

**1. Skip the largest SPAdes tasks in the main run**

Use the `.exitcode` hack with code 100 in the Nextflow work directory to tell Nextflow to treat those tasks as "ignored" (after `errorStrategy` above is configured to ignore code 100).

Example (for each problematic SPAdes work folder):

```bash
echo "100" > <NEXTFLOW_WORK_DIR>/<HASH>/.exitcode
# (repeat for all 9 work dirs)
```

Cancel all running SPAdes single-end jobs so Nextflow can proceed:

```bash
squeue -h -t running -p cpu -u <USERNAME> -o "%i,%.25j" \
  | grep nf-ASSEMBLE_SPADES_SINGLE \
  | cut -f1 -d',' \
  | xargs -i scancel {}
```

**2. Create a separate Nextflow pipeline for these Bioprojects**

```bash
cp -r metatranscriptome metatranscriptome_last9
```

Divide each of the 9 giant Bioprojects into multiple "sub-Bioprojects" for manageability, primarily by:

- Month of submission
- Additional metadata (sample name; single-cell vs spatial) if needed

End result: 74 (later 69) sub-Bioprojects.

**3. Update mapping for last9**

- New `RNaquarium_run_to_bioproject_mapping.tsv` that reflects sub-Bioproject structure.
- Re-create JSON mapping with `create_json_from_mapping.py` in `metatranscriptome_last9`.
- Remove any leftover folders in `all_unmapped_links` that are not in the mapping TSV.
- **Watch out for Windows/MS line breaks:** run `dos2unix` on TSV first.

**4. Include "lost" Bioprojects**

At one resume of the main run, 7 Bioprojects were not recovered correctly and effectively "lost". These 7 were added into the last9 pipeline so that ultimately we end with:

- One main pipeline output.
- One "last9+7" pipeline output.
- Both combined for downstream taxonomic and virus curation steps.

The "last9" pipeline uses fewer chunks (8 nt chunks, 16 Diamond sub-chunks), so per-chunk runtimes are shorter even though SPAdes can be long for the largest sub-bioprojects.

---

## Safely Stopping and Resuming Nextflow on HPC

Graceful stop (recommended):

```bash
# First send INT
scancel --signal=INT <jobid>

# If still running after ~30 seconds
scancel -f --signal=INT <jobid>
```

After fixing configs or code, re-run with `-resume`.

**Caveats:**

Resume can fail if:

- Work directory structure changes (e.g. changed JSON or input folders).
- Container or script paths are modified too aggressively.

We observed resume failures for last9 when folders and JSON no longer matched; that pipeline effectively started from scratch.

---

## Full nt BLAST Performance: `/local/scratch` Strategy

### Problem

The first host-filter BLAST step was OK, but the full nt BLAST (20 chunks) was initially far too slow when run against databases on shared storage:

- Single chunks projected 14–18 days.
- Queue saturation limited throughput despite 50-job Nextflow concurrency.

### Observations

- `core_nt` total directory: ~1.1 TB (including tarballs).
- Actual working BLAST database size: ~250 GB.
- This fits on individual node `/local/scratch`.

### Manual Workaround

**1. Stop main full NT BLASTs within Nextflow**

Cancel 16 running full BLAST chunks (leaving 4 not yet started). This freed up Nextflow's 50-job concurrency for Diamond jobs.

**2. Run full NT BLAST manually with local database**

For each chunk:

- Copy or rsync `core_nt` database to `/local/scratch`.
- Run BLAST with `-db` pointing to `/local/scratch/core_nt`.
- Output written to `/local/scratch` or HPC scratch (local output proved as fast or faster).

Example skeleton for one-off BLAST SLURM scripts:

```bash
DB_SRC="<SHARED_DB_PATH>/2025/core_nt"
DB_DEST="/local/scratch/core_nt"

# Only copy if DB doesn't already exist
if [ ! -d "$DB_DEST/core_nt" ]; then
    echo "Copying core_nt database to node-local disk..."
    mkdir -p "$DB_DEST"
    rsync -a "$DB_SRC"/ "$DB_DEST"/
else
    echo "Database already exists on local disk."
fi

# Run BLAST
<BLASTN_BIN_PATH>/blastn \
  -db "$DB_DEST/core_nt" \
  -query <NEXTFLOW_WORK_DIR>/chunk_004_nonzfhum.fasta \
  -outfmt "6 qseqid sseqid staxid ssciname sskingdom pident length mismatch qcovs gapopen qstart qend sstart send evalue bitscore stitle" \
  -num_threads 4 \
  -evalue 0.05 \
  -max_target_seqs 100 \
  -out /local/scratch/chunk_004_nonzfhum.blast.txt

cp /local/scratch/chunk_004_nonzfhum.blast.txt <NEXTFLOW_WORK_DIR>/chunk_004_nonzfhum.blast.txt
```

### Observed Speed-Up

- Chunks reached ~23% completion in 14 h.
- Extrapolated completion ~2.5–4 days per chunk instead of 14–18 days.
- Effective speed-up greater than 3×, even with fewer CPUs per job (4–8).

### Diamond Performance Note

Based on logs, Diamond NR searches did **not** benefit noticeably from moving DBs to `/local/scratch`. We left Diamond jobs running on shared storage.

---

## Post-Nextflow Processing

Once the Nextflow runs are finished (up to full nt BLAST and Diamond), run the SLURM scripts currently in `src/metatranscriptome/scripts/post-processing/` for post-processing.

These should be moved into your working output folder (the PUB_DIR in run_xxx.sh file), with a subset inside output subfolders:

```
<PUB_DIR>/
<PUB_DIR>/nt_blast
<PUB_DIR>/nr_diamond
```

### Step 0a / 0b: Parse BLASTn and Diamond Outputs with Taxonomizr

Both:

- Parse BLAST / Diamond outputs
- Join against the Taxonomizr database
- Produce comprehensive TSVs with taxonomy fields.

---

## Updated NCBI Taxonomy (Taxonomizr) and NA Handling

### Taxonomizr DB Updates

- Initial runs used a February 2025 Taxonomizr database.
- Some NT hits returned taxids not present in that DB.
- **Fix:** Download updated Taxonomizr database (August 2025) and re-run Steps 0a, 0b, 1, 2, etc.

### New NCBI Taxonomy Ranks and "Superkingdom" Reconstruction

Because of NCBI taxonomy changes (see: [New Ranks in NCBI Taxonomy: Domain & Realm - NCBI Insights](https://ncbiinsights.ncbi.nlm.nih.gov/)):

- "Viruses" now live under the `acellular_root` field.
- Bacteria, Archaea, Eukaryota use the `domain` field.

Taxonomizr's legacy `superkingdom` is effectively reconstructed by:

```r
superkingdom <- dplyr::coalesce(acellular_root, domain)
```

This keeps downstream logic compatible with pre-change taxonomy structure.

### NA Taxonomy Handling

- **Bug:** NA taxonomy entries were previously coerced to "other Eukaryota".
- **Updated behavior:** Any unresolved/NA taxonomy is mapped to `Root_unresolved`.
- **Additional fix:** Rare NA values in `taxname_LCA` are caught and flagged with a note in a `lowcoverage/flag` column.

### Removing Early "Vector" / "Synthetic" Filters

Earlier versions of steps 0a/0b removed entries with "vector" or "synthetic" from `taxname` or `target_title`. This was too aggressive and could drop legitimate taxa (e.g. "Komarekiella delphini-convector"). With improved NA handling and the `Root_unresolved` category, these heuristic filters were removed.

---

## Diamond NR Taxid Recovery via Local NR Database

### Problem

Up to ~10% of Diamond NR hits had `taxid`/`taxname` as `NA` or `"N/A"`:

- Likely due to partial/incomplete NR DB when constructing the Diamond index.
- These hits were formerly dropped entirely from downstream analyses.

### Solution

For each Diamond output file (`*.diamond.txt.gz`):

1. Identify rows where taxid is `NA` or `"N/A"`.
2. Extract unique target IDs.
3. Query local NR BLAST database (`nr_cluster_seq`) via `blastdbcmd` to recover accession, taxid, and species name.

Example:

```bash
blastdbcmd -db nr_cluster_seq -entry "WP_197412162.1" -outfmt "%a %T %S"
# WP_197412162.1 1712383 Rheinheimera sp. EpRS3
```

Scripts in `<PUB_DIR>/nr_diamond/`:

- `slurm_diamond_taxid_update_blastdbcmd.sh` – SLURM wrapper
- `diamond_taxid_update_blastdbcmd.py` – fills missing taxids
- `diamond_taxid_moveold_andrename.py` – moves originals and renames updated files

### Guardrail in Step 0b

Before Diamond taxonomy processing (step 0b), R code checks the first 1000 rows for NA/"N/A" taxids:

```r
if ("taxid" %in% names(single_hit)) {
  taxid_check <- single_hit %>%
    dplyr::slice_head(n = 1000) %>%
    dplyr::pull(taxid)
  has_missing <- any(is.na(taxid_check) | taxid_check == "N/A")
  if (has_missing) {
    stop("WARNING: Your first Diamond dataframe includes some NA values in the taxid column. 
         Please run script slurm_diamond_taxid_update_blastdbcmd.sh in the nr_diamond folder 
         to find taxid values from a local database. Then re-run this script.")
  } else {
    cat("TaxID check passed: no NA or 'N/A' values found in the first 1000 rows.\n")
  }
} else {
  stop("ERROR: 'taxid' column not found in the 'single_hit' dataframe.")
}
```

**Important for reproducibility:**
Before running step 0b, run the Diamond-taxid update script if any NA taxids are detected.

---

## BBDuk Contamination Masking (Post-Search Filter)

### Motivation

Certain NCBI "viral" records are actually adapter/primer contamination, e.g.:

- "Crocidura shantungensis ribovirus 3"
- "Brassica yellows virus"
- Some SARS-CoV-2 records
- Various phage/caudoviricetes entries

These were enriched in:

- The "last9" bioprojects
- Some other bioprojects with heavy internal adapter content

Root cause:

- fastp adapter removal is poor at internal adapters.
- Some contaminated contigs also arise from adapter sequences deposited as "viruses" in NCBI.

Re-running the full metatranscriptome pipeline from scratch with a read-level adapter cleanup would take 3+ weeks. We therefore implemented a post-search BBDuk masking strategy at the contig level.

### Strategy

Rather than re-run the entire metatranscriptome pipeline, we:

1. Applied a post-taxonomy BBDuk mask to the non-host contigs.
2. Identified and excluded contigs with excessive adapter content.
3. Propagated this filter into both:
   - General non-host results
   - Virus-specific outputs and Salmon inputs.

### High-Level Plan

1. Use existing BLAST + Diamond + taxonomy pipeline to generate `taxonomy_hits_nonhost_list.fasta` (all non-host contigs; ~35M contigs, ~15 GB).
2. Run BBDuk to mask adapter/contaminant sequences with `N`.
3. Flag any contig with <150 non-N bases for exclusion.
4. Remove these contigs from all downstream non-host/viral tables via an anti-join.
5. Re-run virus-specific scripts and Salmon quantification on the filtered set.

### BBDuk on Non-Host Contigs (Global Filter)

Starting point: full non-host contigs after initial post-tax steps:

```
taxonomy_hits_nonhost_list.fasta
# ~35M contigs, ~15 GB
```

BBDuk command (run via SLURM script `slurm_2b_bbdukflag.sh`):

```bash
bbduk.sh -Xmx2g threads=auto \
  in=taxonomy_hits_nonhost_list.fasta \
  out=taxonomy_hits_nonhost_list_masked.fasta \
  ref=adapters,artifacts,phix,lambda,pjet,kapa,fastp_adapters_with9added.fasta \
  k=21 mink=10 hdist=1 \
  maskmiddle=f kmask=N fastawrap=0 \
  stats=taxonomy_hits_nonhost_masking_stats.txt
```

Notes:

- `fastp_adapters_with9added.fasta` is a custom adapters file that merges:
  - BBDuk's standard adapters/artifacts
  - fastp adapter set
  - Additional adapters/primers observed in our data (e.g. Clontech primers)

Performance: ~90 seconds for 35M contigs (15 GB) on HPC.

### Flagging Contigs to Exclude (<150bp After Masking)

After masking, remove Ns and flag sequences below 150 bp:

```bash
awk '/^>/ {name=$0; next} {gsub(/N/,""); if(length($0)<150) print name}' \
  taxonomy_hits_nonhost_list_masked.fasta \
  | sed 's/^>//' \
  > taxonomy_hits_nonhost_sequences_to_exclude.txt
```

This text file is then used downstream to filter all joins with contig IDs.

### Integration into R Scripts (Step 1b)

In the R script that combines nt and nr hits (`1_combinecontigs_alluvial_treemap_bbdukfilter.r`):

```r
# Before BBDuk filter
allchunks_diamondnr_andblastntclustered0 <- full_join(
  thresholded_hit_nofishnomammalsNT,
  thresholded_hit_nofishnomammalsNR
)

# BBDuk filter
toremove <- readr::read_tsv("taxonomy_hits_nonhost_sequences_to_exclude.txt",
                            col_names = FALSE) |>
  dplyr::rename(query = X1)

allchunks_diamondnr_andblastntclustered <- dplyr::anti_join(
  allchunks_diamondnr_andblastntclustered0,
  toremove,
  by = "query"
)
```

Virus steps then operate on this BBDuk-filtered non-host set.

---

## Pipeline Sequence Summary

The SLURM scripts are currently located in `src/metatranscriptome/scripts/post-processing/`, and should be moved into your working output folder (the PUB_DIR in run_xxx.sh file), with a subset inside output subfolders.

### With Post-Search BBDuk (Used for 75k Release)

1. **Nextflow:** SPAdes → Host-filter BLAST → Full nt BLAST → Diamond nr
2. **SLURM scripts:** 0a, 0b, 1, 2
3. **BBDuk filtering:** 2b (creates exclusion list), then 1b (re-combine with filter)
4. **Adding sequence & Optional Dark Matter accounting and Metacoder Trees :** 3, 4b, 4t
5. **Virus steps:** v5A (if needed), v5, v6, v7b, v8, v9, v10 (see **Metatranscriptome_virus-steps-Technical-Notes.md**)
6. **Salmon quantification**

### Without BBDuk (see Future Plans)

1. **Nextflow:** SPAdes → Host-filter BLAST → Full nt BLAST → Diamond nr
2. **SLURM scripts:** 0a, 0b, 1, 2, 3, optional 4, 4t
3. **Virus steps:** v5A (if needed), v5, v6, v7, v8, v9, v10
4. **Salmon quantification**

### Future Plans

In a future "clean" full Nextflow run, we plan to integrate BBDuk before SPAdes to handle adapter contamination at the read level rather than post-hoc at the contig level.
