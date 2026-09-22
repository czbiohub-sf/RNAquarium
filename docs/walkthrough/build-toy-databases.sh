#!/bin/bash
# =============================================================================
#  WHAT TO EDIT BEFORE RUNNING — one required, two optional:
#    1. DB_DIR                                      (required) where to build
#    2. VIRUS_NUC_ACC / VIRUS_PROT_ACC / VIRUS_TAXID (optional) the virus to seed
#    3. TAXONOMY_DATE                               (optional) taxonomy vintage
#  Nothing else in this file needs changing.
# =============================================================================
#  build-toy-databases.sh
#
#  *** TEST / POSITIVE-CONTROL ONLY — NOT FOR REAL ANALYSIS ***
#
#  Builds deliberately tiny databases so Part II can be run end to end to validate
#  plumbing, with a declared ground truth. Full explanation + parameter table:
#  part2-test-databases.md. Full databases for real work:
#  https://czbiohub-sf.github.io/RNAquarium/inputs-and-databases.html
#
#  GROUND TRUTH: the only organisms these databases can report are Danio rerio and
#  the one virus you seed. Every other contig is "no hit" BY CONSTRUCTION.
#
#  Everything is built into ONE $DB_DIR because BLAST locates taxdb.btd/.bti by
#  `cd`-ing into the database directory at runtime.
#
#  RUN ON A LOGIN NODE (needs internet), in a build env separate from the pipeline:
#    conda create -n rnaq-dbbuild -c conda-forge -c bioconda blast diamond entrez-direct r-taxonomizr
#    conda activate rnaq-dbbuild
#
#  DIAMOND caveat: the .dmnd format is version-sensitive — build nr with a diamond
#  compatible with the one the pipeline queries with. BLAST v5 is more forgiving.
# =============================================================================
set -eo pipefail

# ---- Where to build the toy databases (EDIT ME) -----------------------------
DB_DIR="/path/to/your/scratch/rnaq-toy-dbs"

# ---- NCBI access ------------------------------------------------------------
# edirect over HTTPS is flaky on many clusters: curl reports "(56) unexpected eof"
# even when the data arrives. `export NCBI_API_KEY=...` before running raises rate
# limits and helps; fetch_or_die below retries and judges success by output produced.
if [ -n "${NCBI_API_KEY:-}" ]; then echo "Using NCBI_API_KEY from environment."; fi

# ---- Host taxa --------------------------------------------------------------
# Tagged at SPECIES taxid; the pipeline's --taxids 7954,9605 are the GENUS ids and
# resolve these by lineage.
DANIO_TAXID="7955"
HUMAN_TAXID="9606"

# ---- The virus we actually want to detect, in BOTH core_nt and nr (optional EDIT) ----
VIRUS_NUC_ACC="OK619588.1"     # nucleotide accession -> toy core_nt
VIRUS_PROT_ACC="WKU61618.1"    # protein accession    -> toy nr
VIRUS_TAXID="3063761"
# Danio rerio nucleotide + protein reps to seed the toy nt/nr (so the 2nd host
# net can fire) are fetched by query below; no need to hardcode accessions.

mkdir -p "$DB_DIR"
cd "$DB_DIR"

# ---- fetch_or_die <db> <query> <format> <out> [cap] -------------------------
# Appends FASTA to <out>. Errors loudly on 0 hits; [cap] limits records to keep the
# toy DBs tiny. Retries, and treats "records written" as success (see curl-56 above).
fetch_or_die() {
    local edb="$1" query="$2" fmt="$3" out="$4" cap="${5:-0}"
    local n
    # hit count is on the UNcapped query (so 0 hits still errors). edirect can emit
    # noise/duplicate lines on a flaky connection, so extract the FIRST pure integer.
    n=$(esearch -db "$edb" -query "$query" 2>/dev/null \
        | xtract -pattern ENTREZ_DIRECT -element Count 2>/dev/null \
        | grep -oE '[0-9]+' | head -1)
    n="${n:-0}"
    if [ "${n:-0}" = "0" ]; then
        echo "ERROR: query returned 0 hits (or the count could not be read) - refusing to build an empty DB." >&2
        echo "       db=$edb  query=[$query]" >&2
        echo "       If you saw curl (56) errors, this may be a flaky NCBI connection, not a bad query:" >&2
        echo "       retry, export NCBI_API_KEY=..., or use a different login node." >&2
        exit 1
    fi
    # normalize cap to a clean integer (0 = no cap)
    cap=$(printf '%s' "$cap" | grep -oE '[0-9]+' | head -1); cap="${cap:-0}"
    if [ "$cap" != "0" ]; then
        echo "  [$edb] '$query' -> $n hit(s) (fetching first $cap)"
    else
        echo "  [$edb] '$query' -> $n hit(s)"
    fi

    local before after tmp attempt
    before=$(grep -c '^>' "$out" 2>/dev/null || true); before="${before:-0}"
    tmp="$(mktemp)"
    for attempt in 1 2 3; do
        : > "$tmp"
        if [ "$cap" != "0" ]; then
            # take the first <cap> UIDs, then fetch just those
            esearch -db "$edb" -query "$query" 2>/dev/null \
                | efetch -format uid 2>/dev/null | grep -oE '[0-9]+' | head -n "$cap" \
                | efetch -db "$edb" -format "$fmt" >> "$tmp" 2>/dev/null || true
        else
            esearch -db "$edb" -query "$query" 2>/dev/null \
                | efetch -format "$fmt" >> "$tmp" 2>/dev/null || true
        fi
        if grep -q '^>' "$tmp" 2>/dev/null; then
            cat "$tmp" >> "$out"; rm -f "$tmp"
            after=$(grep -c '^>' "$out"); after="${after:-0}"
            echo "    fetched $((after - before)) sequence(s) on attempt $attempt"
            return 0
        fi
        echo "    attempt $attempt produced no sequences; retrying after 5s..." >&2
        sleep 5
    done
    rm -f "$tmp"
    echo "ERROR: fetch produced no sequences after 3 attempts (network/NCBI?)." >&2
    echo "       db=$edb query=[$query]. Try: export NCBI_API_KEY=... or a different node." >&2
    exit 1
}

# =============================================================================
# 0. NCBI taxonomy dump (LATEST) + BLAST taxonomy files
#    taxdump feeds DIAMOND makedb (step 3). taxdb.btd/.bti are the universal BLAST
#    lookup — not tied to any DB, just need to sit in $DB_DIR.
#    NOTE: step 4 deliberately overrides nameNode.sqlite with an OLDER pinned dump.
# =============================================================================
wget -N ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir -p taxdump && tar -xzf taxdump.tar.gz -C taxdump      # -> taxdump/names.dmp, taxdump/nodes.dmp (LATEST)
[ -s taxdump/names.dmp ] && [ -s taxdump/nodes.dmp ] \
    || { echo "ERROR: names.dmp/nodes.dmp missing after extracting taxdump.tar.gz" >&2; exit 1; }

wget -N ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar -xzf taxdb.tar.gz                                       # -> taxdb.btd, taxdb.bti
[ -f taxdb.btd ] && [ -f taxdb.bti ] \
    && echo "OK: taxdb.btd + taxdb.bti present in $DB_DIR" \
    || { echo "ERROR: taxdb.btd/.bti missing after extracting taxdb.tar.gz" >&2; exit 1; }

# =============================================================================
# 1. Host-filter DB  (pipeline: --nt_dir $DB_DIR  --nt_db_name host_toy)
#    Danio rerio + Homo sapiens rRNA + mitochondrion ONLY.
#    LIMITATION: catches host by ABUNDANCE (rRNA dominates) but misses most host by
#    DIVERSITY. It shows the step runs; it is not a defensible filter. Expect more
#    host contigs to survive than in a real run — the 2nd net (Danio in the toy
#    nt/nr below) catches them.
# =============================================================================
: > host_toy.fa
: > host_toy.nucl.taxidmap

# Fetch mito (uncapped, ~7) + rRNA (capped to 5) for each species.
# NOTE (verify-on-run): these queries are the REQUIREMENT, not tested commands.
for taxid in "$DANIO_TAXID" "$HUMAN_TAXID"; do
    tmp="host_${taxid}.fa"; : > "$tmp"
    fetch_or_die nuccore "txid${taxid}[Organism:exp] AND mitochondrion[Title] AND refseq[Filter]" fasta "$tmp"
    fetch_or_die nuccore "txid${taxid}[Organism:exp] AND (rRNA[Title] OR \"ribosomal RNA\"[Title]) AND refseq[Filter]" fasta "$tmp" 5
    cat "$tmp" >> host_toy.fa
    # taxid map: every accession fetched for this species -> this species taxid
    grep '^>' "$tmp" | sed 's/^>//; s/ .*//' | awk -v t="$taxid" '{print $1"\t"t}' >> host_toy.nucl.taxidmap
done

makeblastdb -in host_toy.fa -dbtype nucl -parse_seqids \
    -taxid_map host_toy.nucl.taxidmap -out host_toy -title "toy_host_rRNA_mito"
#  -> pipeline: --nt_dir $DB_DIR --nt_db_name host_toy
#     (config default --taxids 7954,9605 resolves 7955/9606 by lineage)

# =============================================================================
# 2. Toy core_nt (nucleotide)  (pipeline: --ntfull_dir $DB_DIR --ntfull_db_name core_nt)
#    Danio rerio reps + the target virus. Danio here lets the downstream 2nd host
#    net (Actinopteri drop) fire AND serves as a fish positive control.
# =============================================================================
: > core_nt_toy.fna
: > core_nt_toy.nucl.taxidmap

# target virus (explicit accession)
efetch -db nuccore -id "$VIRUS_NUC_ACC" -format fasta >> core_nt_toy.fna
printf '%s\t%s\n' "$VIRUS_NUC_ACC" "$VIRUS_TAXID" >> core_nt_toy.nucl.taxidmap

# a few Danio rerio nucleotide RefSeq mRNA reps (capped to 5), guarded
tmp="danio_nuc.fa"; : > "$tmp"
fetch_or_die nuccore "txid${DANIO_TAXID}[Organism:exp] AND biomol_mrna[PROP] AND refseq[Filter]" fasta "$tmp" 5
cat "$tmp" >> core_nt_toy.fna
grep '^>' "$tmp" | sed 's/^>//; s/ .*//' | awk -v t="$DANIO_TAXID" '{print $1"\t"t}' >> core_nt_toy.nucl.taxidmap

makeblastdb -in core_nt_toy.fna -dbtype nucl -parse_seqids \
    -taxid_map core_nt_toy.nucl.taxidmap -out core_nt -title "toy_core_nt"

# =============================================================================
# 3. Toy nr (protein)  (pipeline: --nr_dir $DB_DIR ; db basename "nr" -> nr.dmnd)
#    Danio rerio proteins + the target virus protein. Taxonomy is compiled INTO
#    nr.dmnd at makedb time (--taxonmap/--taxonnodes/--taxonnames); none at runtime.
# =============================================================================
: > nr_toy.faa
: > nr_toy.prot.taxidmap

efetch -db protein -id "$VIRUS_PROT_ACC" -format fasta >> nr_toy.faa
printf '%s\t%s\n' "$VIRUS_PROT_ACC" "$VIRUS_TAXID" >> nr_toy.prot.taxidmap

tmp="danio_prot.fa"; : > "$tmp"
fetch_or_die protein "txid${DANIO_TAXID}[Organism:exp] AND refseq[Filter]" fasta "$tmp" 5
cat "$tmp" >> nr_toy.faa
grep '^>' "$tmp" | sed 's/^>//; s/ .*//' | awk -v t="$DANIO_TAXID" '{print $1"\t"t}' >> nr_toy.prot.taxidmap

# DIAMOND wants an NCBI-style prot.accession2taxid map (accession, accession.version, taxid, gi)
{ printf 'accession\taccession.version\ttaxid\tgi\n'
  awk -F'\t' '{split($1,a,"."); print a[1]"\t"$1"\t"$2"\t0"}' nr_toy.prot.taxidmap
} | gzip > prot.accession2taxid.toy.gz

diamond makedb --in nr_toy.faa -d nr \
    --taxonmap prot.accession2taxid.toy.gz \
    --taxonnodes taxdump/nodes.dmp \
    --taxonnames taxdump/names.dmp
#  -> produces nr.dmnd in $DB_DIR

# =============================================================================
# 4. Taxonomizr nameNode.sqlite  (pipeline: --taxonomy_db $DB_DIR/nameNode.sqlite)
#    *** DELIBERATE: this is the ONLY database built from an OLDER taxonomy dump. ***
#    The public R post-processing bins hits by category using the older taxonomy
#    structure. Built from a CURRENT dump, modern viral-realm lineages
#    (Orthornavirae / Pisuviricota) get binned as "other_Eukaryota" — the virus is
#    still identified correctly, but the viruses table comes out EMPTY.
#    Feb-2025 (the 75k run's vintage) fixes it. Change TAXONOMY_DATE for another.
# =============================================================================
TAXONOMY_DATE="2025-02-01"                 # optional EDIT; pinned to the 75k run's vintage
mkdir -p taxdump_pinned
wget -N "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_${TAXONOMY_DATE}.zip"
unzip -o "taxdmp_${TAXONOMY_DATE}.zip" -d taxdump_pinned    # -> taxdump_pinned/names.dmp, nodes.dmp
[ -s taxdump_pinned/names.dmp ] && [ -s taxdump_pinned/nodes.dmp ] \
    || { echo "ERROR: names.dmp/nodes.dmp missing after unzip taxdmp_${TAXONOMY_DATE}.zip" >&2; exit 1; }

# build nameNode.sqlite from the PINNED dump (overwrite any prior one)
rm -f nameNode.sqlite
Rscript -e '
  suppressMessages(library(taxonomizr))
  read.names.sql("taxdump_pinned/names.dmp", "nameNode.sqlite")
  read.nodes.sql("taxdump_pinned/nodes.dmp", "nameNode.sqlite")
  cat("built nameNode.sqlite from pinned taxonomy dump\n")
'
[ -s nameNode.sqlite ] || { echo "ERROR: nameNode.sqlite was not built" >&2; exit 1; }

echo
echo "===================================================================="
echo " Toy databases built in: $DB_DIR"
echo "   host filter    : host_toy.*        --nt_dir $DB_DIR --nt_db_name host_toy"
echo "   core_nt        : core_nt.*         --ntfull_dir $DB_DIR --ntfull_db_name core_nt"
echo "   nr (diamond)   : nr.dmnd           --nr_dir $DB_DIR   (taxonomy from LATEST dump)"
echo "   taxonomizr     : nameNode.sqlite   --taxonomy_db $DB_DIR/nameNode.sqlite   (pinned ${TAXONOMY_DATE})"
echo "   BLAST taxdb    : taxdb.btd/.bti    (in $DB_DIR; found via cd /db - do NOT pass as params)"
echo
echo " NOTE: nameNode.sqlite is built from the PINNED ${TAXONOMY_DATE} taxonomy dump so the"
echo " public categorization script bins viruses correctly; all other DBs use the latest dump."
echo
echo " GROUND TRUTH (positive control): host = Danio + Homo (rRNA+mito); the ONLY"
echo " reportable non-host organism is virus ${VIRUS_NUC_ACC}/${VIRUS_PROT_ACC} (taxid ${VIRUS_TAXID}),"
echo " plus Danio rerio (2nd host net). EVERYTHING else is dark matter BY DESIGN."
echo " Results here are PLUMBING + POSITIVE CONTROL only, never biology."
echo "===================================================================="
