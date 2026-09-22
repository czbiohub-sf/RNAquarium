#!/bin/bash
#SBATCH --job-name=rnaq-part1-nonhost
#SBATCH --time=21:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=YOUR_PARTITION          # e.g. cpu / general / batch (edit for your site)
#SBATCH -e nonhost_slurm.out/slurm-%j.err
#SBATCH -o nonhost_slurm.out/slurm-%j.out
## Uncomment + set your address to receive email notifications:
##SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --mail-user=your.email@example.org

# =============================================================================
# WHAT TO EDIT BEFORE RUNNING — four things, and nothing else in this file:
#   1. #SBATCH --partition   (top of file)  your cluster's queue name
#   2. PROJECT_DIR           (section 1)    your run area
#   3. conda Case A/B/C      (section 3)    how conda becomes available in a batch job
#   4. --genome-size         (section 5)    your genome's size in bytes
# =============================================================================
# RNAquarium Part I (nonhost) — portable submission template
# -----------------------------------------------------------------------------
# This is a TEMPLATE from docs/walkthrough/. Copy it into your own run area, make
# the four edits above, then SUBMIT FROM THE REPO ROOT pointing at your copy:
#
#     cd /path/to/RNAquarium              # the repo root
#     mkdir -p nonhost_slurm.out          # SLURM needs this to exist BEFORE submit
#     sbatch /path/to/your/run-area/nextflow-submit-part1-nonhost.sh
#
# Nextflow only treats main.nf as a local file when launched from src/nonhost/,
# so this script uses SLURM_SUBMIT_DIR (your dir at submit time = the repo root)
# to find the repo and cd's into src/nonhost itself.
#
# WHAT THIS SCRIPT GUARANTEES
#   * Every run is self-contained: its own timestamped Nextflow log, trace,
#     timeline and report, all written into YOUR RUN AREA (not into the repo).
#   * A fresh run refuses to start on top of a previous run's output, rather
#     than silently mixing results.
#   * The run is VERIFIED when it finishes. Part I can print
#     "Execution complete -- Goodbye" while having produced no non-host reads at
#     all, because the gsnap process resolves to errorStrategy 'ignore' and a
#     failed gsnap is dropped from the channel without failing the workflow.
#     Section 6 turns that silent outcome into an explicit pass/fail.
# =============================================================================

set -eo pipefail    # (deliberately no -u: conda activate can trip 'set -u')

# ---- 1. Paths and options you MUST edit --------------------------------------
# PROJECT_DIR is your RUN AREA — the folder holding your edited params/csv. It is
# NOT the repo (that's found automatically) and NOT the output dir (that's set by
# 'publish-dir:' inside the params YAML). See the directory table in README,
# "Where these files go".
PROJECT_DIR="/path/to/your/run-area"                # e.g. /hpc/scratch/.../rnaq-run
PARAMS_FILE="${PROJECT_DIR}/params.zf-example.yaml" # the copy you edited

# Set RESUME="yes" ONLY to continue an interrupted run that used the SAME work
# dir. Leave it "no" for a first run or any clean re-test. A resumed run is
# allowed to write into an existing output dir; a fresh run is not.
RESUME="no"

# Optional extra config overlay (e.g. a debug overlay). Leave empty for normal use.
EXTRA_CONFIG=""

# ---- 2. Derived paths (no need to edit) --------------------------------------
RUNTAG="$(date +%y%m%d-%H%M%S)"
SCRATCH_DIR="${PROJECT_DIR}/scratch"
WORKDIR="${SCRATCH_DIR}/rnaq-part1-work-zf-example/"
LOGDIR="${PROJECT_DIR}/logs"
# Nextflow's own reportdir default is the RELATIVE string "reports", which lands
# inside the repo at src/nonhost/reports/ and accumulates across runs. Force it
# to an absolute per-run directory in the run area so each run's trace is
# unambiguous and the repo checkout stays clean.
REPORTDIR="${PROJECT_DIR}/reports/part1-${RUNTAG}"
NXFLOG="${LOGDIR}/nextflow_part1_zf-example-${RUNTAG}.log"

# ---- 3. SITE SETUP: make conda available IN THIS BATCH JOB, then activate env -
# A batch job starts a FRESH shell that does NOT inherit the conda you loaded
# interactively, so re-establish it here. This is TWO steps, not one choice:
#
# STEP 3a — make conda available. Do the ONE thing that matches your cluster
#   (the Case A/B/C you identified in README "Prerequisites"):
#     Case A — conda from a cluster module: keep the `module load` line (edit name)
#     Case B — personal Miniforge/Miniconda: comment out `module load`, uncomment `source ~/...`
#     Case C — conda already on PATH: comment out BOTH lines below
module load anaconda                                  # Case A (your site's conda module)
# source ~/miniforge3/etc/profile.d/conda.sh          # Case B (edit path to your install)
#
# STEP 3b — always initialize conda's shell hook, then activate the env from
#   README "Step 1". (Safe/harmless in all three cases; Case A/B still need it so
#   `conda activate` works in a non-interactive shell.)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rnaq-part1                             # env: nextflow=23.10.1, mamba=1.5.*, samtools

# ---- 4. Pre-flight checks ----------------------------------------------------
# Fail here, cheaply, rather than after an hour of index building.
REPODIR="${SLURM_SUBMIT_DIR:-$PWD}"

die() { echo "" >&2; echo "ERROR: $*" >&2; echo "" >&2; exit 1; }

[[ -d "$PROJECT_DIR" ]] || die "PROJECT_DIR does not exist: $PROJECT_DIR"
[[ -f "$PARAMS_FILE"  ]] || die "PARAMS_FILE does not exist: $PARAMS_FILE"
[[ -f "${REPODIR}/src/nonhost/main.nf" ]] || die \
  "main.nf not found under ${REPODIR}/src/nonhost/.
   You must sbatch this script FROM THE REPO ROOT so SLURM_SUBMIT_DIR points there."
command -v nextflow >/dev/null || die "nextflow not on PATH after activating rnaq-part1."

# Read publish-dir and accession-list straight out of the params file, so this
# script and the pipeline can never disagree about where output goes.
PUBLISH_DIR=$(sed -n 's/^publish-dir:[[:space:]]*//p'    "$PARAMS_FILE" | head -1 | sed 's/[[:space:]]*$//')
ACCESSION_CSV=$(sed -n 's/^accession-list:[[:space:]]*//p' "$PARAMS_FILE" | head -1 | sed 's/[[:space:]]*$//')
[[ -n "$PUBLISH_DIR"   ]] || die "could not read 'publish-dir:' from $PARAMS_FILE"
[[ -n "$ACCESSION_CSV" ]] || die "could not read 'accession-list:' from $PARAMS_FILE"
[[ -f "$ACCESSION_CSV" ]] || die "accession-list file does not exist: $ACCESSION_CSV"

# A fresh run must not write on top of a previous run's results.
if [[ "$RESUME" != "yes" && -d "${PUBLISH_DIR}/unmapped_reads" ]]; then
	die "output already exists: ${PUBLISH_DIR}/unmapped_reads
   This is a fresh run (RESUME=\"no\"), so it will not write on top of it.
   Rename the old output aside and resubmit, e.g.:
       mv '${PUBLISH_DIR}' '${PUBLISH_DIR%/}.old-\$(date +%y%m%d-%H%M%S)'
   (Rename rather than delete — a previous run's outputs are useful evidence.)"
fi

mkdir -p "$WORKDIR" "$LOGDIR" "$REPORTDIR" "${REPODIR}/nonhost_slurm.out"
export PATH="$PATH:${REPODIR}/src/nonhost/bin"      # tools from setup-minimal.sh

echo "=== RNAquarium Part I =========================================="
echo "  run tag      : $RUNTAG"
echo "  repo         : $REPODIR"
echo "  params       : $PARAMS_FILE"
echo "  publish dir  : $PUBLISH_DIR"
echo "  work dir     : $WORKDIR"
echo "  nextflow log : $NXFLOG"
echo "  reports      : $REPORTDIR"
echo "  resume       : $RESUME"
echo "================================================================"

# ---- 5. Run ------------------------------------------------------------------
cd "${REPODIR}/src/nonhost"

# Output dir + publish-intermediate come from the params YAML (single source of
# truth). Two flags MUST stay on the command line:
#   --publishQCfiltered : the YAML key 'publish-qc-filtered' silently fails,
#                         because the non-standard 'QC' casing breaks the
#                         kebab-case -> camelCase mapping. Other kebab keys in
#                         the params file (accession-list, publish-dir, ...) map
#                         correctly; this one is a genuine casing edge case.
#   --reportdir         : forces per-run reports into the run area (see §2).
# NOTE ON ARGUMENT ORDER: `-log` and `-c` are LAUNCHER options and must come
# BEFORE the `run` subcommand; `-params-file`, `-profile`, `-work-dir`, `-resume`
# and `--param` are `run` options and come after. This matters: which options each
# subcommand accepts varies by Nextflow version (e.g. `nextflow config` does not
# accept `-params-file` on 23.10.1 but does on later releases). Launcher position
# for `-c` is the form that works across versions.
NXF_ARGS=( -log "$NXFLOG" )
[[ -n "$EXTRA_CONFIG" ]] && NXF_ARGS+=( -c "$EXTRA_CONFIG" )
NXF_ARGS+=(
	run main.nf
	-params-file "$PARAMS_FILE"
	--genome-size 1448808562      # EDIT to your genome's value (awk sum, README Step 6)
	--publishQCfiltered true
	--reportdir "$REPORTDIR"
	-profile slurm,mamba
	-work-dir "$WORKDIR"
	-ansi-log false
)
[[ "$RESUME" == "yes" ]] && NXF_ARGS+=( -resume )

echo "+ nextflow ${NXF_ARGS[*]}"

set +e
NXF_OPTS="-Xms500M -Xmx16G" nextflow "${NXF_ARGS[@]}"
NXF_RC=$?
set -e

# ---- 6. Verify the run actually produced non-host reads ----------------------
# Do not trust the exit status alone. gsnap resolves to errorStrategy 'ignore',
# so a gsnap failure is dropped from the channel and the workflow still exits 0
# with "Execution complete -- Goodbye". Part II's first step then fails with
# "No non-host reads found for <BioProject>", a long way from the real cause.
echo ""
echo "=== Verifying Part I output ===================================="

RUNS=$(awk -F, 'NR==1{for(i=1;i<=NF;i++) if($i=="Run") c=i; next} c && $c!="" {print $c}' "$ACCESSION_CSV")
[[ -n "$RUNS" ]] || die "no 'Run' column values parsed from $ACCESSION_CSV"

MISSING=0
for acc in $RUNS; do
	found=""
	for layout in Paired Single; do
		d="${PUBLISH_DIR}/unmapped_reads/${layout}/${acc}"
		[[ -d "$d" ]] || continue
		# >132 bytes matches the pipeline's own "empty run" threshold
		n=$(find "$d" -maxdepth 1 -name '*.fastq.gz' -size +132c 2>/dev/null | wc -l)
		if [[ "$n" -gt 0 ]]; then
			echo "  OK      $acc  ($n non-host fastq.gz in ${layout}/)"
			found="yes"; break
		fi
	done
	if [[ -z "$found" ]]; then
		echo "  MISSING $acc  (no non-empty *.fastq.gz under unmapped_reads/*/${acc}/)"
		MISSING=$((MISSING+1))
	fi
done

if [[ "$MISSING" -gt 0 || "$NXF_RC" -ne 0 ]]; then
	echo ""
	echo "---------------------------------------------------------------"
	echo "Part I did NOT produce usable non-host reads (nextflow rc=$NXF_RC,"
	echo "$MISSING sample(s) missing). Do NOT start Part II."
	echo ""
	# --reportdir is a pipeline param, and the config builds trace.file from
	# ${params.reportdir} at PARSE time. Whether a command-line --reportdir reaches
	# it is version-dependent, so look in the run area first and then fall back to
	# the launch dir, which is where the config's relative "reports" default lands.
	TRACE=$(ls -t "$REPORTDIR"/trace-*.txt "$PWD"/reports/trace-*.txt 2>/dev/null | head -1)
	if [[ -n "$TRACE" ]]; then
		echo "Trace: $TRACE"
		echo "gsnap task status (status is column 6, exit is column 7):"
		grep -i gsnap "$TRACE" | cut -f1-7 | sed 's/^/    /'
		echo ""
		echo "  status=COMPLETED, exit=0  -> gsnap ran fine; the outputs were"
		echo "     produced but not published. Check the work dir named by the"
		echo "     'hash' column (column 4) -- with publishDir mode 'copy' the"
		echo "     originals are still there -- and grep the log for publish errors:"
		echo "         grep -i publish '$NXFLOG'"
		echo "  status=FAILED/ABORTED     -> gsnap itself failed and was hidden by"
		echo "     errorStrategy 'ignore'. Read .command.log and .command.err in"
		echo "     that work dir."
	else
		echo "No trace file found in $REPORTDIR or $PWD/reports/."
		echo "Search wider with:  find '$REPODIR' '$PROJECT_DIR' -name 'trace-*.txt' -newermt today"
	fi
	echo ""
	echo "Nextflow log: $NXFLOG"
	echo "---------------------------------------------------------------"
	exit 1
fi

echo ""
echo "Part I verified: all $(echo "$RUNS" | wc -w) sample(s) have non-host reads."
echo "Output: ${PUBLISH_DIR}/unmapped_reads/"
echo "Next: build the unmapped_links symlink dir for Part II (README Part II Step 1),"
echo "      and verify the links resolve with 'ls -laL'."
echo "================================================================"
