#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.publishDir = "$PWD"
params.publishIntermediate = true
params.genomeSize = null
params.refIndexes = null
params.cleanupScript = ""
params.tmp = null
params.backupTmp = null
params.backupScratchHack = false
params.nxfUnstageHack = false

params.refGenome = "Danio_rerio.GRCz11.dna_sm.primary_assembly.fa"
params.refGenomeGtf = "Danio_rerio.GRCz11.108.gtf"
params.erccFa = "ERCC92.fa"
params.erccGtf = "ERCC92.gtf"
params.hisatUseTranscript = true

params.retainMixed = true

params.seed = 35854

include {
	hisat2_generate_indexes;
} from './step.0.generate_indexes.nf' params(
	publishDir: params.publishDir,
	hisatUseTranscript: params.hisatUseTranscript,
	refGenome: params.refGenome,
	refGenomeGtf: params.refGenomeGtf,
	erccFa: params.erccFa,
	erccGtf: params.erccGtf,
	tmp: params.tmp,
	backupTmp: params.backupTmp,
	backupScratchHack: params.backupScratchHack,
	nxfUnstageHack: params.nxfUnstageHack
)

params.metaIn = 'step_2_sheet.csv'
params.metaOut = 'step_3_sheet.csv'
include {
	LOAD_METASHEET;
	SAVE_METASHEET;
} from './utils.nf'

process hisat2 {
	label 'hisat2'

	def ALIGNER = "hisat2"
	def SUFFIX_OK = "hisatFiltered.fastq"
	def SUFFIX_NG = "hisatFailed.fastq" // unused for hisat2
	def SAM_NAME = "${ALIGNER}_out.sam"
	def SAM_STAGING = "${ALIGNER}.staging.sam"

	input:
	tuple val(meta), path("m?.fq.gz", arity: '1..2')
	tuple val(idx_basename), path("hisat2_index/*")

	output:
	tuple val(meta), path("Unmapped.out.mate?.*.gz", arity: '1..2'), emit: mates
	//tuple val(meta), path("${SAM_NAME}"), emit: sam
	tuple val(meta), path("${ALIGNER}.stats.txt"), emit: sam_stats
	tuple val(meta), path("stats.txt"), emit: stats
	tuple val(meta), path("metrics.txt"), emit: hisat2_debug

	script:
	// TODO: try lower pen-noncansplice
	// TODO: try -k 10
	// possibly --no-unal --omit-sec-seq
	def ALIGNER_CMD = """hisat2 --seed ${params.seed} \
		--met-file metrics.txt --summary-file stats.txt \
		-p ${task.cpus} -k 5 \
		--pen-noncansplice 12 \
		--sp 1,0 \
		--no-temp-splicesite -t \
		-S ${SAM_STAGING} \
		-x hisat2_index/${idx_basename}"""
	def ALIGNER_CMD_PE = """${ALIGNER_CMD} -1 m1.fq.gz -2 m2.fq.gz"""
	def ALIGNER_CMD_SE = """${ALIGNER_CMD} -U m1.fq.gz"""

	def cond = params.retainMixed ?
		(meta.single_end ? '!flag.secondary && flag.unmap' : '!flag.secondary && (flag.unmap || flag.munmap)') :
		(meta.single_end ? '!flag.secondary && flag.unmap' : '!flag.secondary && (flag.unmap && flag.munmap)')

	def NAMES = "${ALIGNER}_unmapped_names.txt"
	def SAMSTATS_CMD = """samtools view -@ ${task.cpus} ${SAM_STAGING} | cut -f2 | sort | uniq -c > ${ALIGNER}.stats.txt"""
	def GET_NAMES_CMD = """samtools view -@ ${task.cpus} -e '${cond}' ${SAM_STAGING} | cut -f1  > ${NAMES}"""
	def FILTER_CMD = """LC_ALL=C fastq-namefilter ${NAMES} -"""

	"""
	set -euo pipefail
	${!meta.single_end ? ALIGNER_CMD_PE : ALIGNER_CMD_SE}

	${SAMSTATS_CMD}
	${GET_NAMES_CMD}
	for i in ${!meta.single_end ? "{1..2}" : "1"} ; do
		${task.ext.gzipCmd} -kcd m\${i}.fq.gz | ${FILTER_CMD} | ${task.ext.gzipCmd} -nc > Unmapped.out.mate\${i}.${SUFFIX_OK}.gz.staging
	done
	mv ${SAM_STAGING} ${SAM_NAME}
	for file in *.gz.staging ; do mv \$file \${file/.gz.staging/.gz} ; done

	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
}

def ensure_hisat2_indexes(ref_indexes,
						  ref_genome, ref_genome_gtf, ercc, ercc_gtf) {
	if (ref_indexes
		&& (indexes = file(ref_indexes))
		&& indexes.exists()) {
		return [indexes.listFiles()[0].getSimpleName(), file(indexes.resolve("*.{ht2,ht2l}")) ]
	} else {
		indexes = hisat2_generate_indexes(file(ref_genome),
										  file(ref_genome_gtf),
										  file(ercc),
										  file(ercc_gtf))
	}
	return indexes
}


workflow {
	indexes = ensure_hisat2_indexes(params.refIndexes,
									params.refGenome,
									params.refGenomeGtf,
									params.ercc,
									params.erccGtf)

	fastqs = LOAD_METASHEET(params.metaIn)

	unmapped_reads = hisat2(fastqs, indexes)

	SAVE_METASHEET(unmapped_reads, params.metaOut)
}
