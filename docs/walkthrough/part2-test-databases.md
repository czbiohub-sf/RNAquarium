# Part II Test Databases — a smoke-test shortcut

> ⚠️ **TEST / SMOKE-TEST ONLY — NOT FOR REAL ANALYSIS.** ⚠️
>
> The databases described here are deliberately tiny. They exist so you can run
> the Part II metatranscriptome pipeline **end to end once** and confirm the
> *plumbing* works — assembly, the BLAST/DIAMOND searches, taxonomy parsing, and
> taxonomy assignment all run, and the output directory structure is correct.
>
> Results are meaningful **only for the handful of sequences you deliberately put
> in** (here, one target virus). Everything else is reported as "no hit / dark
> matter". For any real analysis, build the full databases — see
> [Inputs & Database Requirements](https://czbiohub-sf.github.io/RNAquarium/inputs-and-databases.html).

---

## Why this is possible

Part II runs standard tools — `blastn` against the nt / core_nt databases,
`diamond blastx` against NR, and `taxonomizr` for taxid → lineage lookups. None of
these care whether a database is 1 TB or 1 MB; they just need a database in the
right *format* at the path you give them. So we build the same formats with the
same tools (`makeblastdb`, `diamond makedb`) from a handful of sequences instead
of all of GenBank.

## About the built-in `test` profile

`src/metatranscriptome/nextflow.config` **does** define a `test` profile — but it
only sets database parameter *defaults to placeholder paths*
(`/path/to/databases/...`), along with `process.cache = 'lenient'`,
`num_blast_chunks = 4`, and `num_diamond_chunks = 5`. It does **not** ship small
databases. So the toy databases here are still needed. You can either:

- pass the toy database paths directly on the command line, or
- run `-profile test,slurm,singularity` and override just the DB params with the toy
  paths on the command line (CLI overrides the test-profile defaults).

**Which dependency profile does Part II use?** Part II runs its tools in **containers
via `-profile slurm,singularity`** — that's the proven Biohub setup, and it's what
`run-part2-metatranscriptome.sh` uses. With `-profile singularity`, Nextflow pulls
and runs the containers for you; you never invoke `singularity` yourself, but it must
be available on your cluster (usually `module load singularity`/`apptainer`). The
config *also* defines `conda` and `mamba` profiles, but those haven't been tested for
Part II here — stick with singularity unless you have a reason not to.

The `test` profile also *declares* three taxonomy params — `tax_btd`, `tax_bti`,
`tax4blast` — but **no process actually reads them** on public main, and the real
`run_75k.sh` does not pass them. They're vestigial. BLAST instead finds its taxonomy
files (`taxdb.btd`/`taxdb.bti`) because they sit **inside the database directory** that
gets mounted at `/db`; the BLAST step `cd`s there at runtime. So the requirement is not
"pass a param" but "put `taxdb.btd`/`taxdb.bti` in the same folder as the db" — which the
builder does. `taxonomy4blast.sqlite3` is **not** needed.

## The real databases, for reference

The full databases are large (from the Inputs & Databases doc):

| Database | Approx. size | How the docs build it |
| -------- | ------------ | --------------------- |
| core_nt (full BLASTn) | ~1.1 TB | `update_blastdb.pl --decompress core_nt` |
| NR (DIAMOND) | ~300 GB | `wget nr.gz` → `diamond makedb --in nr.gz -d nr` |
| Taxonomizr SQLite | ~75 GB* | `prepareDatabase("nameNode.sqlite")` in R |
| BLAST taxdb (`.btd`/`.bti`) | < 1 GB | from `taxdb.tar.gz`, placed in each db dir |
| nt host-only subset | ~250 GB | clustered subset of nt |

\* The ~75 GB is dominated by the accession→taxid table. The pipeline's
`--taxonomy_db` uses only the much smaller **names/nodes** (`nameNode.sqlite`)
part, which is what the toy builder produces in full (so taxonomy still resolves
correctly for whatever taxids your toy hits carry).

## Tools you need to BUILD the toy databases

These are **not** the pipeline's tools — the pipeline supplies its own per step via
`-profile mamba` (Part I) / `-profile singularity` (Part II), so you do **not** add them to the
`rnaq-part1` / `rnaq-part2` Nextflow envs. They're only needed for this one-time,
by-hand build step: `makeblastdb` / `update_blastdb.pl` (blast), `diamond`,
`efetch` (entrez-direct), and `Rscript` + `taxonomizr`.

Easiest is a dedicated build env with everything from conda, including R + `taxonomizr`:

```bash
conda create -n rnaq-dbbuild -c conda-forge -c bioconda blast diamond entrez-direct r-taxonomizr
conda activate rnaq-dbbuild
```

The builder needs `makeblastdb`/`diamond`/`efetch` on `PATH` and an `Rscript` that can
`library(taxonomizr)`. If conda's R misbehaves on your cluster, substitute your own R
module with `taxonomizr` installed (see the walkthrough README's Known Issues).

> **Version caveat:** build the DIAMOND database with a diamond version
> **compatible with the one the pipeline uses to query it** — the DIAMOND db
> format is version-sensitive, and a db built by a newer diamond can be unreadable
> by an older one. BLAST v5 databases are more forgiving. If a search step
> complains about database format, match your build tool's version to the
> pipeline's (check the metatranscriptome env/container after cloning).

## Building the toy databases

Run `build-toy-databases.sh` (in this folder) on a **login node** (it needs
internet). Edit the paths and the virus accessions at the top first:

- `DB_DIR` — where to build (all toy databases go in this one folder)
- `VIRUS_NUC_ACC` / `VIRUS_PROT_ACC` / `VIRUS_TAXID` — the target virus to make
  detectable (defaults: `OK619588.1` / `WKU61618.1` / taxid `3063761`)

The Danio + Homo host sequences are fetched by NCBI query (guarded to error if a query
returns nothing), so you don't hardcode accessions.

It produces, in `DB_DIR`:

| Toy database | Files | Pipeline parameter |
| ------------ | ----- | ------------------ |
| Host filter | `host_toy.*` | `--nt_dir $DB_DIR --nt_db_name host_toy` |
| core_nt | `core_nt.*` | `--ntfull_dir $DB_DIR --ntfull_db_name core_nt` |
| NR (DIAMOND) | `nr.dmnd` | `--nr_dir $DB_DIR` |
| Taxonomizr | `nameNode.sqlite` | `--taxonomy_db $DB_DIR/nameNode.sqlite` |
| BLAST taxdb | `taxdb.btd` / `taxdb.bti` | *(not a param — must sit in `$DB_DIR`)* |

How each toy database is made, briefly:

- **Host filter (`host_toy`)** — Danio rerio + Homo sapiens **rRNA + mitochondrion**
  only, tagged with taxids `7955` / `9606`. The pipeline's `--taxids 7954,9605` (genus
  Danio / Homo) resolves these by lineage, so no override is needed. This is a *minimal*
  host DB (rRNA+mito), enough to demonstrate the host-filter step on the toy data — not a
  substitute for the full host database a real run uses.
- **core_nt / NR** — Danio rerio reps **plus** the target virus, each built with taxids
  attached (`makeblastdb -taxid_map` for core_nt; DIAMOND's
  `--taxonmap/--taxonnodes/--taxonnames` for NR) so column-4 scientific names resolve —
  **without this, the R post-processing silently drops every row** (see the walkthrough
  README's "Part II: taxonomy tables come back empty").
- **Taxonomizr** — built from `names.dmp` + `nodes.dmp` only (the small part), not the
  giant accession→taxid table. The builder downloads the **latest** taxonomy dump for
  everything, **except** it builds `nameNode.sqlite` (the file the categorization step
  reads) from a pinned **Feb-2025** dump — a single, deliberate override. A current dump
  makes the public taxonomy script miscategorize viral-realm lineages
  (`Orthornavirae`/`Pisuviricota`) as `other_Eukaryota` instead of `viruses`; the Feb-2025
  vintage (what the public docs say the 75k run used) fixes it and makes the run
  reproducible. Change `TAXONOMY_DATE` in the script for a different vintage.
- **BLAST taxdb** — the real `taxdb.btd`/`taxdb.bti` from `taxdb.tar.gz`, placed in
  `$DB_DIR` so BLAST finds them via `cd /db`. Not passed as a param.

## What a toy run does and does not validate

**Validates (the point):** the pipeline assembles the Part I non-host reads, runs
each search step without error, parses taxonomy through taxonomizr (PROCESS_BLAST /
PROCESS_DIAMOND) and draws the alluvial plot, and writes the expected output structure.
(Combining nt+nr, curation, and Salmon quantification are post-Nextflow SLURM/R scripts,
out of scope here.) All on a machine without
1.7 TB of databases. If your target virus is present in the reads, you should see
it recovered and labeled with taxid `3063761`, which is a nice sanity check that
detection works end to end.

**Does not validate:** anything about organisms you did not put in the toy
databases. Do not read taxonomy calls, LCA assignments, or abundances for anything
other than your seeded virus. Swap in the full databases before drawing any real
conclusion.
