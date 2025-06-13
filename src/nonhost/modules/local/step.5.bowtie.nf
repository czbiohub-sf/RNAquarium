#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.publishDir = "$PWD"
params.publishIntermediate = true
params.genomeSize = null
params.bowtieIndex = null
params.cleanupScript = ""
params.tmp = null
params.backupTmp = null
params.backupScratchHack = false
params.nxfUnstageHack = false

params.refGenome = "Danio_rerio.GRCz11.dna_sm.primary_assembly.fa"
params.refGenomeGtf = "Danio_rerio.GRCz11.108.gtf"
params.erccFa = "ERCC92.fa"
params.erccGtf = "ERCC92.gtf"

params.retainMixed = true

params.seed = 32854

include {
	bowtie2_generate_indexes;
} from './step.0.generate_indexes.nf' params(
	publishDir: params.publishDir,
	refGenome: params.refGenome,
	refGenomeGtf: params.refGenomeGtf,
	erccFa: params.erccFa,
	erccGtf: params.erccGtf,
	tmp: params.tmp,
	backupTmp: params.backupTmp,
	backupScratchHack: params.backupScratchHack,
	nxfUnstageHack: params.nxfUnstageHack
)

params.metaIn = 'step_4_sheet.csv'
params.metaOut = 'step_5_sheet.csv'
include {
	LOAD_METASHEET;
	SAVE_METASHEET;
} from './utils.nf'

process bowtie2 {
	label 'bowtie2'

	def ALIGNER = "bowtie2"
	def SUFFIX_OK = "filteredbyBT.fastq"
	def SUFFIX_NG = "failedBT.fastq" // unused for bowtie2
	def SAM_NAME = "${ALIGNER}_out.sam"

	input:
	tuple val(meta), path("m?.fq.gz", arity: '1..2')
	tuple val(idx_basename), path("bowtie2_index/*")

	output:
	tuple val(meta), path("Unmapped.out.mate?.*.gz", arity: '1..2'), emit: mates
	//tuple val(meta), path("${SAM_NAME}"), emit: sam
	tuple val(meta), path("${ALIGNER}.stats.txt"), emit: stats

	script:
	//def index_base = file(index_dir).listFiles()[0].getSimpleName()
	def ALIGNER_CMD = """bowtie2 --quiet --very-sensitive-local -p $task.cpus \
		--rg-id na --rg LB:na --rg SM:na --rg PL:na --rg PU:na \
		--seed ${params.seed} -x bowtie2_index/${idx_basename} """
	def ALIGNER_CMD_PE = """${ALIGNER_CMD} -1 m1.fq.gz -2 m2.fq.gz -S ${ALIGNER}.staging.sam"""
	def ALIGNER_CMD_SE = """${ALIGNER_CMD} -U m1.fq.gz -S ${ALIGNER}.staging.sam"""

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
	${!meta.single_end ? ALIGNER_CMD_PE : ALIGNER_CMD_SE}
	mv ${ALIGNER}.staging.sam ${SAM_NAME}

	${SAMSTATS_CMD}
	${GET_NAMES_CMD}
	for i in ${!meta.single_end ? "{1..2}" : "1"} ; do
		${task.ext.gzipCmd} -kcd m\${i}.fq.gz | ${FILTER_CMD} | ${task.ext.gzipCmd} -nc > Unmapped.out.mate\${i}.${SUFFIX_OK}.gz.staging
	done
	for file in *.gz.staging ; do mv \$file \${file/.gz.staging/.gz} ; done

	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
}

def ensure_bowtie2_indexes(ref_indexes,
						   ref_genome, ercc) {
	if (ref_indexes
		&& (indexes = file(ref_indexes))
		&& indexes.exists()) {
		return [indexes.listFiles()[0].getSimpleName(), file(indexes.resolve("*.{bt2,bt2l}"))]
	} else {
		indexes = bowtie2_generate_indexes(ref_genome,
										   file(ercc))
	}
	return indexes
}
