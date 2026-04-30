# Metatranscriptome Pipeline (Part II)

Taxonomic classification and viral discovery from non-host reads.

## Features

* Input: non-host reads from Part I (transcriptomic + filtering pipeline)
* Assembles contigs with SPAdes
* BLASTn against nt and Diamond against nr for taxonomic classification
* BBDuk-based adapter/contaminant filtering
* Virus-specific curation with realm taxonomy support
* Salmon quantification of viral contigs
* Distributed computing on HPC environments (Slurm job scheduler)

## Pipeline Architecture

**Important:** The Nextflow pipeline only handles assembly and database searches. All downstream processing is performed by standalone SLURM/R scripts.

| Phase | Tool | Steps |
|-------|------|-------|
| Assembly & Search | Nextflow | SPAdes assembly → Host-filter BLAST → Full nt BLAST → Diamond nr |
| Taxonomy & Filtering | SLURM scripts | Taxonomy parsing, contig filtering, BBDuk masking |
| Virus Curation | SLURM scripts | Realm taxonomy, clustering, metacoder trees |
| Quantification | SLURM scripts | Salmon viral contig quantification |

All post-Nextflow scripts are located in `src/metatranscriptome/scripts/post_processing/`.

## Usage

### Clone the repository

```bash
git clone https://github.com/czbiohub-sf/Zebrafish-RNAquarium
```

### Install dependencies

The metatranscriptome pipeline requires the following dependencies in `PATH`:

* SPAdes ( __https://github.com/ablab/spades__ )
* BLAST+ ( __https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/__ )
* Diamond ( __https://github.com/bbuchfink/diamond__ )
* BBTools/BBDuk ( __https://jgi.doe.gov/data-and-tools/software-tools/bbtools/__ )
* Salmon ( __https://combine-lab.github.io/salmon/__ )
* SeqKit ( __https://bioinf.shenwei.me/seqkit/__ )
* samtools ( __https://www.htslib.org/download/__ )
* Taxonomizr R package ( __https://github.com/sherrillmix/taxonomizr__ )

Alternatively, run the Nextflow pipeline with `-profile conda`, `-profile mamba`, `-profile docker`, or `-profile singularity` to automatically pull containers or environments for most steps.

### Prepare inputs

The metatranscriptome pipeline consumes unmapped reads from Part I. Ensure all unmapped reads (single-end and paired-end) are symlinked into a single directory:

```bash
cd Zebrafish-RNAquarium/src/metatranscriptome
mkdir -p all_unmapped_links
# Symlink your unmapped reads here
```

Create a run-to-bioproject mapping file (`RNaquarium_run_to_bioproject_mapping.tsv`) and generate the JSON mapping:

```bash
python scripts/create_json_from_mapping.py
```

### Database requirements

* **NCBI nt database** (core_nt) for BLASTn searches
* **NCBI nr database** for Diamond searches
* **Taxonomizr database** (recommend August 2025 or later for realm support)

### Running the pipeline

**Phase 1: Nextflow (assembly & searches)**

```bash
sbatch run_75k.sh
```

Or with Nextflow directly:

```bash
nextflow run main.nf \
  --input_dir all_unmapped_links \
  --mapping_json bioproject_mapping.json \
  --nt_db /path/to/core_nt \
  --nr_db /path/to/nr \
  -profile slurm
```

**Phase 2: SLURM scripts (taxonomy & downstream)**

After Nextflow completes, run the post-processing scripts sequentially. See the technical notes for the complete pipeline sequence.

## Documentation

* **docs/Metatranscriptome-Technical-Notes.md** – full pipeline documentation, workarounds, and step-by-step SLURM script sequence
* **docs/Metatranscriptome_virus-steps-Technical-Notes.md** – virus-specific curation and Salmon quantification

## Outputs

| Output | Description |
|--------|-------------|
| Non-host contigs | SPAdes-assembled contigs after host filtering |
| Taxonomy tables | BLASTn (nt) and Diamond (nr) results with full taxonomy |
| BBDuk-filtered contigs | Adapter-masked contigs with contamination removed (post-search BBDuk, used for 75k release) |
| Virus tables | Curated viral contigs with realm taxonomy |
| Salmon counts | Quantification matrix (~75k runs × ~190k viral contigs) |

## Key considerations

* **Large bioprojects**: May need to be subdivided; see technical notes for "last9" strategy
* **BLAST performance**: Place nt database on `/local/scratch` for 3×+ speedup
* **Taxonomy updates**: Use Taxonomizr DB from August 2025+ for proper viral realm support
* **Diamond NA taxids**: Run `slurm_diamond_taxid_update_blastdbcmd.sh` if NA taxids detected

See **docs/Metatranscriptome-Technical-Notes.md** for complete details.
