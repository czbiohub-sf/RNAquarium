# Virus-Specific Curation & Salmon Quantification

## Virus LCA and Taxonomy

Virus-specific steps operate on non-host contigs classified as viral by either nt or nr (LCA logic).

Key points:

- NCBI changed viral top-level rank from "clade" to "realm"; we now explicitly add realm in virus scripts.
- Virus scripts rely on:
  - `taxname_lca_NTorNR`
  - `taxid_lca_NTorNR` (added for portal-friendliness)
  - Realm and other ranks reconstructed via Taxonomizr from updated DB.

### Virus SLURM/R Pipeline Sequence

**With Post-Search BBDuk (used for 75k release):**

(v5A if needed) → v5 → v6 → **v7b** → v8 → v9 → v10

**Without BBDuk (original sequence):**

(v5A if needed) → v5 → v6 → v7 → v8 → v9 → v10

| Step | Description |
|------|-------------|
| v5A | Add realm info to virus taxonomy – **only required if realm data is missing** (see guardrail below) |
| v5 | Initial virus contig filtering and annotation |
| v6 | Clustering and deduplication |
| v7 | Metacoder tree generation |
| v7b | Metacoder tree generation **with BBDuk filter** – filters & trims viral sequences before clustering (use instead of v7 when using BBDuk filter & Step 1b) |
| v8 | Non-vertebrate / allbutchordates subset creation |
| v9 | Final virus table assembly |
| v10 | Salmon quantification prep |

### Guardrail in Step v5: Viral Realm Taxonomy Check

NCBI changed the viral top-level rank from "clade" to "realm", so older Taxonomizr databases may have all-NA values in the `tax_clade` columns for viruses. Step v5 checks whether viral realm names (Riboviria, Varidnaviria, Duplodnaviria) are present in the clade columns:

```r
# Define the expected viral realm names
viral_realms <- c("Riboviria", "Varidnaviria", "Duplodnaviria")

# Check if any of these realm names exist in either clade column
realms_in_ntclustered <- any(allchunks_diamondnr_andblastntclustered$tax_clade_NTclustered %in% viral_realms, na.rm = TRUE)
realms_in_nr <- any(allchunks_diamondnr_andblastntclustered$tax_clade_NR %in% viral_realms, na.rm = TRUE)

if (!realms_in_ntclustered && !realms_in_nr) {
  stop("WARNING: Your virus dataframe is missing the highest level taxonomy, 
       please run script slurm_virus_5A_curation.sh to pull the new 'realm' field 
       to replace your missing values in 'clade'.")
} else {
  cat("Viral realm taxonomy found - proceeding with analysis.\n")
}
```

**If the check fails:** Run `slurm_virus_5A_curation.sh` (step v5A) to pull the new "realm" field using Taxonomizr, then re-run step v5.

**If the check passes:** Step v5A can be skipped; proceed directly with v5.

**Outputs produced:**

- Virus contig tables
- Clusters
- Metacoder trees
- Subsets for non-vertebrate / allbutchordates etc.

Additional general metacoder trees (fungi, etc.) are produced by:

- SLURM: `slurm_4t_general_metacoder_trees.sh`
- R script: `general_metacoder_trees.r`

(Outputs saved under `RNQuarium_outputs/metacoder_trees`.)

---

## BBDuk on Virus Contigs for Salmon

For Salmon quantification, we use a virus-only FASTA with adapters masked and removed, plus minimum length enforcement.

### Starting Points

Virus sequences and clusters (no nt targets), from 09-17 analysis

```
taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-11.fasta
```

Also removing all Sprivivirus contigs and replacing with 2 reference genomes:

```
taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-17_sprivireplacediwth2genomes.fasta
```

### BBDuk Mask and Hard Removal of Ns

```bash
# On HPC
module load anaconda
conda activate bbmap

bbduk.sh -Xmx2g threads=auto \
  in=taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-17_sprivireplacediwth2genomes.fasta \
  out=taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-17_sprivireplacediwth2genomes_masked.fa \
  ref=adapters,artifacts,phix,lambda,pjet,kapa,fastp_adapters_with9added.fasta \
  k=21 mink=10 hdist=1 \
  maskmiddle=f kmask=N fastawrap=0 \
  stats=taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-17_masking_stats.txt
```

Then remove Ns:

```bash
sed '/^>/!s/N//g' \
  taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-17_sprivireplacediwth2genomes_masked.fa \
  > taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-17_sprivireplacediwth2genomes_masked_trimmed0.fa
```

### Re-ordering, Phage/Non-Phage Grouping, Min Length

Using SeqKit:

```bash
module load anaconda
conda activate SeqKit

# Split phage/prokaryote vs non-phage
seqkit grep -r -i -p "phage|Caudovi|prokaryote" -w 0 \
  taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-17_sprivireplacediwth2genomes_masked_trimmed0.fa \
  > phage.fa

seqkit grep -v -r -i -p "phage|Caudovi|prokaryote" -w 0 \
  taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-17_sprivireplacediwth2genomes_masked_trimmed0.fa \
  > nonphage.tmp.fa

# Separate "special" contaminated groups
seqkit grep -r -i -p "respiratory_syndrome|Crocidura|carnation|porcine" -w 0 \
  nonphage.tmp.fa > nonphage_special.fa

seqkit grep -v -r -i -p "respiratory_syndrome|Crocidura|carnation|porcine" -w 0 \
  nonphage.tmp.fa > nonphage_regular.fa

# Sort each group
seqkit sort -n -w 0 nonphage_regular.fa > nonphage_regular.sorted.fa
seqkit sort -n -w 0 nonphage_special.fa > nonphage_special.sorted.fa
seqkit sort -n -w 0 phage.fa > phage.sorted.fa

# Concatenate in desired order
cat nonphage_regular.sorted.fa nonphage_special.sorted.fa phage.sorted.fa \
  > taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-17_sprivireplacediwth2genomes_masked_trimmed1.fa

# Enforce min length 150 bp
seqkit seq -m 150 -w 0 \
  taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-17_sprivireplacediwth2genomes_masked_trimmed1.fa \
  > reordered_trimmed_min150.fa

# Final cleanup of Ns (should be none)
sed '/^>/!s/N//g' reordered_trimmed_min150.fa \
  > taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-17_sprivireplacediwth2genomes_masked_trimmed.fa
```

This final FASTA is what we feed into Salmon:

```
<OUTPUT_ROOT>/salmon_steps/nonhost_viruscounts_salmon/taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-17_sprivireplacediwth2genomes_masked_trimmed.fa
```

Companion TSV (without sequence, no nt targets) for annotations, with updates using h5ad:

```
<OUTPUT_ROOT>/salmon_steps/nonhost_viruscounts_salmon/taxonomy_hits_viruses_withclusters_nosequencenotargetsforsalmon_sprivireplacediwth2genomes_fishassociated_2025-09-17.tsv
<OUTPUT_ROOT>/adata_obj/75k_unstable_anndata_zfin_aliases_metadata.h5ad
```

---

## Salmon Commands

Working directory:

```
<OUTPUT_ROOT>/host_mapping/unmapped_reads/nonhost_viruscounts_salmon
```

We used the following SLURM scripts:

- `slurm_salmonpreSEPE_sept24.sh`
- `slurm_salmonPE_sept24.sh`
- `slurm_salmonSE_sept24.sh`
- `slurm_salmonSEb_sept24.sh` (split SE into two batches)
- `slurm_salmonSEPEpost_sept24.sh` (post-processing & aggregation)

Final Salmon count matrix (~75k SRA runs × ~182k virus contigs):

```
<OUTPUT_ROOT>/host_mapping/unmapped_reads/all.quantt_sept24.tsv.gz
```

This is the file passed to downstream ML work.

---

## Data Portal-Specific Fields

For downstream data portal use, we added:

- `taxid_lca_NTorNR`: a single, easily parseable LCA taxid derived from nt/nr LCAs.
- Non-vertebrate and allbutchordates subsets from the main taxonomy tables.

These changes are integrated into Steps 1 and 1b R scripts.

---

## Future Plan (Clean Nextflow Run)

For a future "clean" metatranscriptome run (e.g. for publication-grade reruns), we plan to:

1. **Integrate BBDuk directly into the Nextflow pipeline:**
   - At the very start of the metatranscriptome pipeline (before SPAdes)
   - Only masking/removing adapter regions
   - Not at the host-counting stage (to avoid inflating runtime by 39,000× reads)

2. **Use the updated August 2025 Taxonomizr database** from the beginning.

3. **Run with an updated Diamond nr database** if available, ideally one that doesn't exhibit the large NA taxid issue.

4. **Keep the `/local/scratch` nt BLAST strategy** as the default.
