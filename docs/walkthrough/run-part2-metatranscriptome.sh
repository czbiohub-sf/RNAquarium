#!/usr/bin/env bash
#SBATCH --job-name=rnaq-part2-metatx
#SBATCH --time=21:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=YOUR_PARTITION          # e.g. cpu / general / batch (edit for your site)
#
# EDIT THESE THREE PATHS to your repo (SLURM #SBATCH lines can't use variables, so they
# must be literal). --chdir starts the job IN the pipeline dir so `main.nf` and the
# relative container path './containers/taxonomy.sif' resolve. -o/-e are ABSOLUTE so
# SLURM can write them regardless of where you submit from (a relative path here is
# resolved at submit time, before the job runs — an easy footgun).
#SBATCH --chdir=/path/to/RNAquarium/src/metatranscriptome
#SBATCH -o /path/to/RNAquarium/metatranscriptome_slurm.out/slurm-%j.out
#SBATCH -e /path/to/RNAquarium/metatranscriptome_slurm.out/slurm-%j.err
## Uncomment + set your address to receive email notifications:
##SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --mail-user=your.email@example.org

# =============================================================================
# WHAT TO EDIT BEFORE RUNNING:
#   1. #SBATCH --partition          (top of file)  your cluster's queue name
#   2. #SBATCH --chdir / -o / -e    (top of file)  three LITERAL paths (no variables)
#   3. PROJECT_DIR                  (section 1)    your run area
#   4. DB_DIR                       (section 2)    your databases (toy or full)
#   5. conda Case A/B/C             (section 3)    how conda becomes available in a job
# =============================================================================
# RNAquarium Part II (metatranscriptome) — portable submission template
# -----------------------------------------------------------------------------
# This is a TEMPLATE from docs/walkthrough/. Copy it into your run area, make the
# five edits above, then create the SLURM log dir once and submit:
#     mkdir -p /path/to/RNAquarium/metatranscriptome_slurm.out   # must exist before submit
#     sbatch /path/to/your/run-area/run-part2-metatranscriptome.sh
#
# REQUIRES, before you run this:
#   1. Part I finished — its unmapped_reads/Paired/ outputs exist.
#   2. The --unmerged_accessions links dir built from those outputs (see README).
#   3. Reference databases — full ones, OR the toy DBs from build-toy-databases.sh
#      (see part2-test-databases.md). Toy DBs validate PLUMBING ONLY, not biology.
#   4. The taxonomy.sif container built + the test-scale conf/modules.config swapped in
#      (README Part II Step 4).
#
# --chdir (above) makes the job start in src/metatranscriptome/, so Nextflow finds
# main.nf and the './containers/taxonomy.sif' path used by the taxonomy steps.
#
# NOTE: Part II parameter names are snake_case and the output flag is --publish_dir
#       (underscore) — different from Part I's kebab-case --publish-dir.
# =============================================================================
set -eo pipefail

# ---- 1. Paths you MUST edit -------------------------------------------------
PROJECT_DIR="/path/to/your/run-area"                # e.g. /hpc/scratch/.../rnaq-run
SCRATCH_DIR="${PROJECT_DIR}/scratch"                # Nextflow work dir; on scratch, outside the repo
BIOPROJ_MAP="${PROJECT_DIR}/zf-example-bioproject-mapping.json"
UNMERGED_ACC="${PROJECT_DIR}/unmapped_links/"       # dir of symlinks you built (see README)
PUB_DIR="${PROJECT_DIR}/RNAquarium_outputs/metatranscriptome_zf-example/"
WORK_DIR="${SCRATCH_DIR}/rnaq-part2-work-zf-example/"
LOGDIR="${PROJECT_DIR}/logs/"

# ---- 2. Reference databases (edit to your site) -----------------------------
# For the toy test DBs from build-toy-databases.sh, set DB_DIR to that folder and
# the names below already match what the builder produces. These are the ONLY database
# params the public metatranscriptome pipeline reads (matching run_75k.sh).
DB_DIR="/path/to/your/databases"                    # e.g. .../rnaq-toy-dbs for a toy run
NT_DIR="${DB_DIR}";        NT_DB_NAME="host_toy"     # host-filter db (full run: nt_compressed_shuffled.fa)
NTFULL_DIR="${DB_DIR}";    NTFULL_DB_NAME="core_nt"  # full BLASTn db
NR_DIR="${DB_DIR}"                                   # DIAMOND nr (basename nr -> nr.dmnd)
TAXONOMY_DB="${DB_DIR}/nameNode.sqlite"             # taxonomizr names/nodes
# NOTE: tax_btd/tax_bti/tax4blast and taxonomy_db_fallback are NOT used by the public
# pipeline (they only appear in its `test` profile) — run_75k.sh does not pass them, so
# neither do we. The toy builder still makes taxdb.btd/.bti; they're just not passed here.

# ---- 3. SITE SETUP: make conda + singularity available IN THE JOB -----------
# A batch job starts a FRESH shell that does NOT inherit your interactive modules.
# conda is needed to run Nextflow (from the rnaq-part2 env); Singularity is needed
# because Part II runs its tools in containers.
#
# Conda is TWO steps, not one choice (same as Part I):
#   2a) make conda available — the ONE thing matching your cluster (Case A/B/C):
module load anaconda                                  # Case A (your site's conda module)
# source ~/miniforge3/etc/profile.d/conda.sh          # Case B (edit path); Case C: comment both
#   2b) always init the hook + activate:
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rnaq-part2                             # env: nextflow=24.10.5

# Container runtime (Apptainer/Singularity) — Part II runs its tools in containers,
# so a container runtime must be on PATH IN THE JOB. This varies by cluster:
#   * Check yours once, interactively:  command -v apptainer ; command -v singularity
#   * If EITHER prints a path, it's already on PATH — you need NO module line here
#     (leave the line below commented). This is the common case.
#   * If NEITHER prints, module-load it. The module name varies (apptainer,
#     singularity, or a versioned name); find it with `module avail apptainer` or
#     `module spider apptainer`, then uncomment and edit the line below:
# module load apptainer        # <- UNCOMMENT + set to your cluster's name, ONLY if needed
#
# (The `nextflow run` line below uses `-profile slurm,singularity`. Apptainer usually
# provides a `singularity` command alias, so that profile works either way. If your
# cluster has apptainer but NO `singularity` command, change it to `slurm,apptainer`.)
#
# Part II uses containers (-profile singularity below). Nextflow pulls and runs the
# containers for you — you never invoke singularity yourself. This matches the
# proven Biohub Part II setup. (conda/mamba profiles also exist but are untested
# here; singularity is the recommended, verified path.)
unset NXF_VER                                        # ignore any stale pinned version

# ---- 4. Make output dirs ----------------------------------------------------
# --chdir already put us in src/metatranscriptome/ (so main.nf and
# ./containers/taxonomy.sif resolve). Just create the output/work/log dirs.
mkdir -p "${WORK_DIR}" "${PUB_DIR}" "${LOGDIR}"

# Isolated cache so Part II can run at the same time as a Part I job
export NXF_CACHE_DIR="${WORK_DIR}/.nextflow_cache"

REPORT="${PUB_DIR}/report.zf-example.html"
TIMELINE="${PUB_DIR}/timeline.zf-example.html"

# ---- 5. Run (already in src/metatranscriptome/ via --chdir) -----------------
nextflow -log "${LOGDIR}/nextflow_part2_zf-example.log" run \
    -profile slurm,singularity \
    -resume \
    -with-report "${REPORT}" \
    -with-timeline "${TIMELINE}" \
    -w "${WORK_DIR}" \
    -ansi-log false \
    --bioproj_map "${BIOPROJ_MAP}" \
    --unmerged_accessions "${UNMERGED_ACC}" \
    --publish_dir "${PUB_DIR}" \
    --nt_dir "${NT_DIR}" \
    --nt_db_name "${NT_DB_NAME}" \
    --ntfull_dir "${NTFULL_DIR}" \
    --ntfull_db_name "${NTFULL_DB_NAME}" \
    --nr_dir "${NR_DIR}" \
    --taxonomy_db "${TAXONOMY_DB}" \
    --num_blast_chunks 2 \
    --num_diamond_chunks 2 \
    main.nf
