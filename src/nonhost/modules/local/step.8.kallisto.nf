#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.publishDir = "$PWD"
params.publishIntermediate = true
params.genomeSize = null
//params.kbContamIndexes = null
params.cleanupScript = ""
params.tmp = null
params.backupTmp = null
params.backupScratchHack = false
params.nxfUnstageHack = false

params.kbRetainMixed = false

params.seed = 32854

include {
	kb_generate_indexes;
} from './step.0.generate_indexes.nf' params(
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
)

process kb_negative {
	label 'kb'

	def ALIGNER = "kallisto"
	def SUFFIX_OK = "kallistoFiltered.fastq"
	def SUFFIX_NG = "kallistoFailed.fastq"
	//def BUS_NAME = "${ALIGNER}_out.bus"
	//def BUS_STAGING = "${ALIGNER}.staging.bus"

	input:
	tuple val(meta), path("m?.fq.gz", arity: '1..2')
	path(kb_contam_index)

	output:
	tuple val(meta), path("Unmapped.out.mate?.*.gz", arity: '1..2'), emit: mates
	tuple val(meta), path("run_info.json"), emit: stats

	// for discard-mixed: default behavior
	// for keep-mixed: must run independently on paired, intersect nums
	// note: kallisto returns an error (1) exit code on 0 reads aligned.
	script:
	def ALIGNER_CMD = """kallisto bus -i "${kb_contam_index}" \
		-x bulk -t ${task.cpus} --num -o . """
	def ALIGNER_CMD_PE = """${ALIGNER_CMD} --paired m1.fq.gz m2.fq.gz"""
	def ALIGNER_CMD_SE = """${ALIGNER_CMD} m1.fq.gz"""
	def FILTER_CMD = "fastq-numfilter nums.idx \$total_reads -"

	if (meta.single_end || !params.kbRetainMixed)
	"""
	set -euo pipefail
	# kallisto exit code is 1 when 0 reads aligned, too.
	${!meta.single_end ? ALIGNER_CMD_PE : ALIGNER_CMD_SE} || [ \$? -eq 1 ]

	# grab number of input reads
	total_reads=\$(grep -m1 "n_processed" run_info.json | grep -om1 "[0-9]\\+")
	# convert output.bus to spot indices
	bustools text -f -p output.bus | cut -f5 > nums.idx

	for i in ${!meta.single_end ? "{1..2}" : "1"} ; do
		${task.ext.gzipCmd} -kcd m\${i}.fq.gz | ${FILTER_CMD} | ${task.ext.gzipCmd} -nc > Unmapped.out.mate\${i}.${SUFFIX_OK}.gz.staging
	done
	for file in *.gz.staging ; do mv \$file \${file/.gz.staging/.gz} ; done

	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
	else
	"""
	set -euo pipefail
	# kallisto exit code is 1 when 0 reads aligned, too.
	${ALIGNER_CMD} m1.fq.gz || [ \$? -eq 1 ]
	bustools text -f -p output.bus | cut -f5 > nums1.idx	
	${ALIGNER_CMD} m2.fq.gz || [ \$? -eq 1 ]
	bustools text -f -p output.bus | cut -f5 | sort nums1.idx - | uniq -d > nums.idx
	total_reads=\$(grep -m1 "n_processed" run_info.json | grep -om1 "[0-9]\\+")

	for i in ${!meta.single_end ? "{1..2}" : "1"} ; do
		${task.ext.gzipCmd} -kcd m\${i}.fq.gz | ${FILTER_CMD} | ${task.ext.gzipCmd} -nc > Unmapped.out.mate\${i}.${SUFFIX_OK}.gz.staging
	done
	for file in *.gz.staging ; do mv \$file \${file/.gz.staging/.gz} ; done

	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
}


def ensure_kb_indexes(index,
					  fa) {
	if (index
		&& (index_file = file(index))
		&& index_file.exists()) {
		return index_file
	} else if (fa) {
		indexes = kb_generate_indexes(files(fa))
		return indexes
	}
	return null
}
