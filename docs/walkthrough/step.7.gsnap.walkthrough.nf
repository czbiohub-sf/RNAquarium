#!/usr/bin/env nextflow
// ============================================================================
// WALKTHROUGH SWAP 5 — derived from src/nonhost/modules/local/step.7.gsnap.nf
//
// Copy over the repo file (git-revertible, like the other walkthrough swaps):
//     cp step.7.gsnap.walkthrough.nf \
//        <repo>/src/nonhost/modules/local/step.7.gsnap.nf
//
// ONE change, in the script block. Everything else is byte-identical to
// upstream main. See the "[walkthrough swap 5]" comment below.
// ============================================================================

nextflow.enable.dsl=2

params.publishDir = "$PWD"
params.publishIntermediate = true
params.genomeSize = null
params.maxMismatch = 0.3
params.cleanupScript = ""
params.tmp = null
params.backupTmp = null
params.backupScratchHack = false
params.nxfUnstageHack = false

params.metaIn = 'step_6_sheet.csv'
params.metaOut = 'step_7_sheet.csv'
include {
	LOAD_METASHEET;
	SAVE_METASHEET;
} from './utils.nf'

include {
	gsnap_generate_indexes;
} from './step.0.generate_indexes.nf' params(
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
)

process gsnap {
	label 'gmap'

	def ALIGNER = "gsnap"
	def SUFFIX_OK = "filteredbyBT.dedup.gsnapFiltered.fastq"
	def SUFFIX_NG = "filteredbyBT.dedup.gsnapSkipped.fastq"
	def SAM_NAME = "${ALIGNER}_out.sam"

	input:
	tuple val(meta), path("m?.fq.gz", arity: '1..2')
	path index_dir

	output:
	tuple val(meta), path("Unmapped.out.mate?.*.gz", arity: '1..2'), emit: mates
	tuple val(meta), path("${SAM_NAME}"), emit: sam
	tuple val(meta), path("${ALIGNER}.stats.txt"), emit: stats

	script:
	def gsnap_gmap_bin = params.genomeSize < 2**32 ? "gsnap" : "gmap"
	def index_name = file(index_dir).getName()
	def ALIGNER_CMD = """${gsnap_gmap_bin} -A sam \
		--batch=4 \
		--use-shared-memory=0 \
		--maxsearch=128 \
		--npaths=1 `# maximum paths to print` \
		--ordered -t ${task.cpus} \
		--max-mismatches ${params.maxMismatch} \
		-D . -d ${index_name} \
		-o ${ALIGNER}.staging.sam """
	def ALIGNER_CMD_PE = """${ALIGNER_CMD} m1.fq m2.fq"""
	def ALIGNER_CMD_SE = """${ALIGNER_CMD} m1.fq"""

	// filter settings
	def PRIMARY = '!flag.secondary && !flag.supplementary'
	def cond = params.retainMixed ?
		(meta.single_end ? "${PRIMARY} && flag.unmap" : "${PRIMARY} && (flag.unmap || flag.munmap)") :
		(meta.single_end ? "${PRIMARY} && flag.unmap" : "${PRIMARY} && (flag.unmap && flag.munmap)")
	def NAMES = "${ALIGNER}_unmapped_names.txt"
	def SAMSTATS_CMD = """samtools view -@ ${task.cpus} ${SAM_NAME} | cut -f2 | sort | uniq -c > ${ALIGNER}.stats.txt"""
	def GET_NAMES_CMD = """samtools view -@ ${task.cpus} -e '${cond}' ${SAM_NAME} | cut -f1  > ${NAMES}"""
	def FILTER_CMD = """LC_ALL=C fastq-namefilter ${NAMES} -"""

	"""
	${!meta.single_end
	? "${task.ext.gzipCmd} -kcd m1.fq.gz > m1.fq ; ${task.ext.gzipCmd} -kcd m2.fq.gz > m2.fq"
	: "${task.ext.gzipCmd} -kcd m1.fq.gz > m1.fq"
	}
	set +e  # suppress terminate-on-error
	${!meta.single_end ? ALIGNER_CMD_PE : ALIGNER_CMD_SE}
	gsnap_rc=\$?  # [walkthrough swap 5] capture IMMEDIATELY -- see note below
	set -e  # resume terminate on error, check error and clear outfile.

	# =====================================================================
	# [walkthrough swap 5] Upstream reads:
	#
	#     if [[ \$? > 0 ]]; then
	#
	# That test can never be true, for two independent reasons:
	#   1. \$? is the exit status of the preceding `set -e`, not of gsnap --
	#      gsnap's status was discarded one line earlier. (The May version of
	#      this file captured `gsnap_rc=\$?` immediately, as restored above.)
	#   2. Inside [[ ]], `>` is a STRING comparison, not numeric, so even a
	#      correct status would compare lexicographically.
	#
	# The result is that gsnap can never be detected as having failed, and the
	# gsnapSkipped fallback below -- whose entire purpose is to pass the dedup
	# reads through unchanged when gsnap fails -- is unreachable. A failed
	# gsnap therefore yields empty outputs, a task that still exits 0, and a
	# pipeline that reports "Execution complete -- Goodbye" having silently
	# discarded every sample. Part II then fails far downstream with
	# "No non-host reads found for <BioProject>".
	#
	# The `-s` test is the second half of the fix. gsnap has been observed
	# exiting 0 while producing an empty or absent SAM (and, given empty input,
	# it exits 0 and writes no -o file at all). Treating "no usable SAM" as
	# failure catches that case regardless of exit status, so the fallback
	# fires and the run still produces usable reads.
	# =====================================================================
	if [[ \$gsnap_rc -gt 0 || ! -s ${ALIGNER}.staging.sam ]]; then  # gsnap failed, pass through reads unchanged...
		echo "gsnap FAILED (rc=\$gsnap_rc, staging sam \$(stat -c%s ${ALIGNER}.staging.sam 2>/dev/null || echo missing) bytes) -- passing reads through unfiltered" >&2
		:> ${SAM_NAME}
		:> ${ALIGNER}.stats.txt
		for i in ${!meta.single_end ? "{1..2}" : "1"} ; do
			cp m\${i}.fq.gz Unmapped.out.mate\${i}.${SUFFIX_NG}.gz.staging
		done
	else
		mv ${ALIGNER}.staging.sam ${SAM_NAME}
		${SAMSTATS_CMD}
		${GET_NAMES_CMD}
		for i in ${!meta.single_end ? "{1..2}" : "1"} ; do
			<m\${i}.fq ${FILTER_CMD} | ${task.ext.gzipCmd} -nc > Unmapped.out.mate\${i}.${SUFFIX_OK}.gz.staging
		done
	fi
	for file in *.gz.staging ; do mv \$file \${file/.gz.staging/.gz} ; done

	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
}

def ensure_gsnap_indexes(ref_indexes,
						 ref_genome, ercc) {
	if (ref_indexes
		&& (indexes = file(ref_indexes))
		&& indexes.exists()) {
	} else {
		indexes = gsnap_generate_indexes(ref_genome,
										 file(ercc))
	}
	return indexes
}
