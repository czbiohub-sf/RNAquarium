# Part II Test Databases — a smoke-test shortcut

> ⚠️ **TEST ONLY — NOT FOR REAL ANALYSIS.** These databases are deliberately tiny.
> They let you run Part II end to end once and confirm the *plumbing* works.
> For real work build the full databases —
> [Inputs & Database Requirements](https://czbiohub-sf.github.io/RNAquarium/inputs-and-databases.html).

Part II's tools (`blastn`, `diamond blastx`, `taxonomizr`) don't care whether a
database is 1 TB or 1 MB — only that it's in the right format at the path you give it.
So the builder makes the same formats with the same tools from a handful of sequences.

---

## Build them

```bash
conda create -n rnaq-dbbuild -c conda-forge -c bioconda blast diamond entrez-direct r-taxonomizr
conda activate rnaq-dbbuild

# edit the top of build-toy-databases.sh first:
#   DB_DIR                                        — where to build (one folder for all)
#   VIRUS_NUC_ACC / VIRUS_PROT_ACC / VIRUS_TAXID  — the virus to make detectable
#                                                   (default OK619588.1 / WKU61618.1 / 3063761)
bash build-toy-databases.sh        # LOGIN NODE — needs internet
```

Host sequences (Danio + Homo) are fetched by NCBI query, so no accessions to hardcode.

If conda's R misbehaves, use your own R with `taxonomizr` instead (see the walkthrough
README's Known Issues).

> **DIAMOND version:** build the NR database with a diamond version compatible with the
> one that will query it — the `.dmnd` format is version-sensitive. BLAST v5 is more
> forgiving.

## What you get

Everything lands in one `DB_DIR`, because BLAST finds its taxonomy files by `cd`-ing
into the mounted database directory.

| Toy database | Files | Pipeline parameter |
| ------------ | ----- | ------------------ |
| Host filter | `host_toy.*` | `--nt_dir $DB_DIR --nt_db_name host_toy` |
| core_nt | `core_nt.*` | `--ntfull_dir $DB_DIR --ntfull_db_name core_nt` |
| NR (DIAMOND) | `nr.dmnd` | `--nr_dir $DB_DIR` |
| Taxonomizr | `nameNode.sqlite` | `--taxonomy_db $DB_DIR/nameNode.sqlite` |
| BLAST taxdb | `taxdb.btd` / `taxdb.bti` | *not a param — must sit in `$DB_DIR`* |

Do **not** pass `--tax_btd` / `--tax_bti` / `--tax4blast`. They are declared but no
process reads them, and `run_75k.sh` doesn't use them. `taxonomy4blast.sqlite3` isn't
needed either.

## Two things the builder does deliberately

**Taxids are attached to every database** (`makeblastdb -taxid_map` for core_nt,
DIAMOND's `--taxonmap/--taxonnodes/--taxonnames` for NR). Without them, column 4 comes
back `N/A` and **the R post-processing silently drops every row** — you get empty
taxonomy tables and no error. See the README's "Part II: taxonomy tables come back
empty".

**`nameNode.sqlite` is built from a pinned Feb-2025 taxonomy dump**, while everything
else uses current NCBI data. With a current dump the public taxonomy script
miscategorizes viral-realm lineages (`Orthornavirae`/`Pisuviricota`) as
`other_Eukaryota` instead of `viruses`. Feb-2025 — the vintage the 75k run used — fixes
that and makes the run reproducible. Change `TAXONOMY_DATE` in the script for a
different vintage.

## What a toy run proves

**It validates:** assembly runs, every search step completes, taxonomy parses through
taxonomizr, the alluvial plot draws, and the output structure is right — without 1.7 TB
of databases. As a positive control, a single virus should come back labelled with
its taxid (~3 virus contigs in the example).

**It does not validate:** anything about organisms you didn't seed. Ignore taxonomy
calls, LCA assignments and abundances for everything else — those are artifacts of a
tiny database, not biology.

---

## Notes

**The built-in `test` profile doesn't ship databases.** `src/metatranscriptome/nextflow.config`
defines a `test` profile, but it only sets DB parameters to placeholder paths plus
`cache = 'lenient'` and chunk counts — so you still need these toy databases. Either
pass the toy paths on the command line, or use `-profile test,slurm,singularity` and
override the DB params (CLI wins).

**Use `-profile slurm,singularity` for Part II.** Nextflow pulls and runs the containers
itself; you never invoke singularity directly, but it must exist on the cluster. The
`conda`/`mamba` profiles also exist but are untested for Part II.

**For scale, the real databases:**

| Database | Approx. size |
| -------- | ------------ |
| core_nt (BLASTn) | ~1.1 TB |
| NR (DIAMOND) | ~300 GB |
| Taxonomizr SQLite | ~75 GB |
| nt host-only subset | ~250 GB |
| BLAST taxdb | < 1 GB |

The 75 GB is mostly the accession→taxid table; `--taxonomy_db` uses only the much
smaller names/nodes portion, which the toy builder produces in full — so taxonomy
resolves correctly for whatever taxids your hits carry.
