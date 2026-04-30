# Inputs & Database Requirements

This page consolidates all input files, reference data, and external databases needed to run the RNAquarium pipeline. Part I (Transcriptomic + Filtering) and Part II (Metatranscriptomics) have distinct requirements — Part II consumes the unmapped reads produced by Part I.

---

## Part I: Transcriptomic + Filtering

### Host reference genome (required)

A primary genome assembly for the host organism. Must be **bgzip-compressed** and **indexed with `samtools faidx`**.

For zebrafish, the 75k production run used **GRCz12tu** (NCBI RefSeq):

| File | Example |
|------|---------|
| Genome FASTA | `GCF_049306965.1_GRCz12tu_genomic.fa.gz` |
| GTF annotation | `GCF_049306965.1_GRCz12tu_genomic.gtf` |
| Genome size (bp) | `1448808683` |

Download from NCBI: <https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_049306965.1/>

The previous release used Ensembl GRCz11 release 108 (`Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz`). Although the code defaults still reference GRCz11, the 75k run overrode these via a params file — see [Parameters](parameters.md) for details.

Preparing the genome:

```bash
# If starting from an uncompressed FASTA:
bgzip genome.fa
samtools faidx genome.fa.gz
```

### GTF annotation file (required)

A standard GTF file matching the genome assembly. Used for gene-level read counting (STAR/featureCounts/HTSeq).

For GRCz12tu: `GCF_049306965.1_GRCz12tu_genomic.gtf` (from the same NCBI Datasets link above).

### ERCC spike-in sequences (required)

External RNA Controls Consortium (ERCC) synthetic spike-in sequences, used for QC normalization.

| File | Description |
|------|-------------|
| `ERCC92.fa` | ERCC spike-in FASTA sequences |
| `ERCC92.gtf` | ERCC annotation in GTF format |

Download from Thermo Fisher: <https://tools.thermofisher.com/content/sfs/manuals/ERCC92.zip>

### SRA accession list (required)

A text file with one SRA run accession per line (e.g., `SRR1234567`). The pipeline will download FASTQ data directly from NCBI SRA.

Example:

```
SRR1234567
SRR1234568
SRR1234569
```

A test list is provided at `src/nonhost/data/SRA_accession_list.test.txt`. The 75k production run used a list of 77,188 zebrafish RNA-seq accessions.

### Contaminant sequences (optional)

A FASTA file containing known contaminant sequences (e.g., common plasmid vectors, RefSeq contaminants). Used by the kallisto contamination filter (Step 2.5).

The 75k run used `drerio_refseq_and_addgene_regions.fa`, which combines Danio rerio RefSeq sequences with Addgene plasmid regions.

Related parameter: `--contam-fa`

### Aligner indexes (optional, auto-generated)

If not provided, the pipeline generates indexes automatically from the genome FASTA. Pre-built indexes save significant time on repeated runs. Five sets are needed:

| Aligner | Parameter | Example directory name |
|---------|-----------|----------------------|
| HISAT2 | `--hisat-ref-indexes` | `hisat2_GCF_049306965.1_GRCz12tu_genomic_indexes/` |
| STAR | `--star-ref-indexes` | `star_GCF_049306965_indexes/` |
| STAR + ERCC | `--star-ref-indexes-ercc` | `star_GCF_049306965_indexes.ERCC/` |
| Bowtie2 | `--bowtie-ref-indexes` | `bowtie2_GCF_049306965.1_GRCz12tu_genomic_indexes/` |
| GMAP/GSNAP | `--gsnap-ref-indexes` | `gmap_GCF_049306965_indexes/` |

Optionally, a kallisto contamination index: `--kb-contam-indexes`

---

## Part II: Metatranscriptomics

### Unmapped reads from Part I (required)

The primary input to Part II is the set of non-host reads produced by Part I. These are gzipped FASTQ files organized into `Single/` and `Paired/` subdirectories under `unmapped_reads/`.

Parameter: `--unmerged_accessions` (path to the directory containing symlinks to all unmapped read files)

### Run-to-bioproject mapping (required)

A JSON file that maps SRA run accessions to their parent BioProjects. This is used to organize assembly into logical groups.

Generated from NCBI metadata using `scripts/create_json_from_mapping.py`, which reads the SRA Accessions metadata table:

| File | Description |
|------|-------------|
| `SRA_Accessions.tab` | NCBI SRA metadata dump |
| `bioproject_mapping.json` | Generated JSON mapping (pipeline input) |

The SRA metadata table can be downloaded from: `ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab`

### NCBI NT database — host-only subset (required)

Used for the initial host-filtering BLAST step. A clustered/compressed subset of NT is sufficient for this step since it only checks for host taxids.

| Parameter | Description |
|-----------|-------------|
| `--nt_dir` | Directory containing the BLAST database |
| `--nt_db_name` | Database name (e.g., `nt_compressed_shuffled.fa`) |

The 75k run used a clustered NT subset. Host taxids are defined by the `--taxids` parameter (default: `7954,9605` for Danio rerio and human).

BLAST parameters for this step: evalue `1e-20`, max_target_seqs `1`.

### NCBI NT database — full core_nt (required)

Used for the comprehensive BLASTn search against all assembled contigs that pass the host filter.

| Parameter | Description |
|-----------|-------------|
| `--ntfull_dir` | Directory containing the full core_nt BLAST database |
| `--ntfull_db_name` | Database name (e.g., `core_nt`) |

Download from NCBI:

```bash
# Download the core_nt database (approximately 250 GB compressed, ~1.1 TB total)
mkdir core_nt && cd core_nt
update_blastdb.pl --decompress core_nt
```

BLAST parameters for this step: evalue `0.05`, max_target_seqs `100`.

**Performance note:** Placing the NT database on local scratch storage (`/local/scratch`) provides ~3x speedup over network-mounted filesystems.

### NCBI NR database (required)

Used for Diamond blastx searches for protein-level taxonomic classification.

| Parameter | Description |
|-----------|-------------|
| `--nr_dir` | Directory containing the Diamond-formatted NR database |

Download and prepare:

```bash
# Download NR
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

# Build Diamond database
diamond makedb --in nr.gz -d nr
```

Diamond search parameters: `blastx --ultra-sensitive --top 3` (with `qcovhsp` in output format).

### Taxonomizr database (required)

An SQLite database used by the R `taxonomizr` package for NCBI taxonomy lookups. This maps taxids to full lineage information (superkingdom, phylum, class, order, family, genus, species) and newer ranks (domain, realm, clade).

| Parameter | Description |
|-----------|-------------|
| `--taxonomy_db` | Path to `nameNode.sqlite` |

Build from within R:

```r
library(taxonomizr)
prepareDatabase("nameNode.sqlite")
```

**Version recommendation:** Use a database built from **August 2025 or later** NCBI taxonomy dumps for proper viral realm support. The initial 75k Nextflow run used a February 2025 database; post-processing scripts were subsequently re-run with the August 2025 database to resolve missing taxids and incorporate updated viral taxonomy (see [Metatranscriptome Technical Notes](Metatranscriptome-Technical-Notes.md) for details).

### BLAST taxonomy files (required for host-filter step)

Additional NCBI taxonomy files used by BLAST for taxid-based filtering:

| Parameter | File | Description |
|-----------|------|-------------|
| `tax_btd` | `taxdb.btd` | BLAST taxonomy database (binary) |
| `tax_bti` | `taxdb.bti` | BLAST taxonomy index |

Download from: `ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz`

### BBDuk adapter/contaminant sequences (used in post-processing)

A merged FASTA file of adapter and contaminant sequences used by BBDuk for post-search contamination masking. The file `fastp_adapters_with9added.fasta` combines standard adapter sets (Illumina, Clontech) with custom additions.

Located at: `src/metatranscriptome/fastp_adapters_with9added.fasta` (included in the repository)

---

## Database Summary Table

| Database | Pipeline Part | Approx. Size | Update Frequency | Required |
|----------|--------------|-------------|------------------|----------|
| Host genome + GTF | Part I | ~1.5 GB | Per assembly release | Yes |
| ERCC spike-ins | Part I | < 1 MB | Rarely | Yes |
| NT (clustered subset) | Part II | ~250 GB | Quarterly | Yes |
| NT (full core_nt) | Part II | ~1.1 TB | Quarterly | Yes |
| NR (Diamond-formatted) | Part II | ~300 GB | Quarterly | Yes |
| Taxonomizr (SQLite) | Part II | ~75 GB | As needed | Yes |
| BLAST taxdb files | Part II | < 1 GB | With NT/NR updates | Yes |
| Contaminant FASTA | Part I | < 10 MB | Rarely | Optional |
| BBDuk adapters | Part II | < 1 MB | Rarely | Included |

---

## Disk Space Planning

A full RNAquarium deployment requires substantial storage:

- **Databases:** ~1.7 TB minimum (NT + NR + Taxonomizr + taxonomy files)
- **Aligner indexes (Part I):** ~50-100 GB for all five mapper indexes
- **Working space (Part I):** Depends on number of SRA runs; the 75k run required multiple TB of scratch space
- **Working space (Part II):** Assembly and BLAST results scale with the number of non-host contigs; the 75k run produced outputs in the TB range
- **Local scratch:** Strongly recommended for BLAST databases to improve I/O performance
