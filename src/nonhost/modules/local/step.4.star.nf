#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.publishDir = "$PWD"
params.publishIntermediate = true
params.genomeSize = null
params.refIndexesErcc = null // "Danio_rerio.GRCz11.108.ERCC"
params.refIndexes = null //"Danio_rerio.GRCz11.108"
params.cleanupScript = ""
params.tmp = null
params.backupTmp = null
params.backupScratchHack = false
params.nxfUnstageHack = false

params.refGenome = "Danio_rerio.GRCz11.dna_sm.primary_assembly.fa"
params.refGenomeGtf = "Danio_rerio.GRCz11.108.gtf"
params.erccFa = "ERCC92.fa"
params.erccGtf = "ERCC92.gtf"

params.starUseSharedMem = false
params.starThreadsSmall = 4
params.starThreadsLarge = 16

params.retainMixed = true

params.seed = 32854

include {
	star_generate_indexes;
} from './step.0.generate_indexes.nf' params(
	publishDir: params.publishDir,
)

params.metaIn = 'step_3_sheet.csv'
params.metaOut = 'step_4_sheet.csv'
include {
	LOAD_METASHEET;
	SAVE_METASHEET;
} from './utils.nf'


process star {
	label 'star'

	def ALIGNER = "star"
	def SUFFIX_OK = "STAR.fastq"
	def SUFFIX_NG = "failedSTAR.fastq" // unused for STAR
	def SAM_NAME = "${ALIGNER}_out.bam"

	input:
	tuple val(meta), path("m?.fq.gz", arity: '1..2')
	path indexes_dir // "Danio_rerio.GRCz11.108.ERCC"

	
	output:
	tuple val(meta), path("Unmapped.out.mate?.*.gz", arity: '1..2'), emit: mates
	//tuple val(meta), path("${SAM_NAME}"), emit: sam
	tuple val(meta), path("Log.final.out"), emit: stats
	tuple val(meta), path("${ALIGNER}.stats.txt"), emit: sam_stats

	script:
	// note: we could get read num with --outSAMreadID Number
	//	--seedSearchStartLmax 30 \
	//	--alignSJoverhangMin 3 \
	def loadType = params.starUseSharedMem ? "LoadAndRemove" : "NoSharedMemory"
	def ALIGNER_CMD = """STAR --outFilterMultimapNmax 99999 \
		--limitOutSJcollapsed 200000000 \
		--outFilterMismatchNmax 999 \
		--outFilterScoreMinOverLread 0.5 \
		--outFilterMatchNminOverLread 0.5 \
		--outFilterType Normal \
		--outSAMtype BAM Unsorted \
		--outSAMmode NoQS \
		--outSAMunmapped Within `# needed for samstats` \
		--outSAMmultNmax 1 \
		--runRNGseed ${params.seed} \
		--genomeLoad $loadType \
		--genomeDir ${indexes_dir} \
		--readFilesCommand ${task.ext.gzipCmd} -dc \
		--runThreadN ${task.cpus} """
	def ALIGNER_CMD_PE = """${ALIGNER_CMD} --readFilesIn m1.fq.gz m2.fq.gz"""
	def ALIGNER_CMD_SE = """${ALIGNER_CMD} --readFilesIn m1.fq.gz"""

	// filter settings
	def cond = params.retainMixed ?
		(meta.single_end ? '!flag.secondary && flag.unmap' : '!flag.secondary && (flag.unmap || flag.munmap)') :
		(meta.single_end ? '!flag.secondary && flag.unmap' : '!flag.secondary && (flag.unmap && flag.munmap)')

	def NAMES = "${ALIGNER}_unmapped_names.txt"
	def SAMSTATS_CMD = """samtools view -@ ${task.cpus} ${SAM_NAME} | cut -f2 | sort | uniq -c > ${ALIGNER}.stats.txt"""
	def GET_NAMES_CMD = """samtools view -@ ${task.cpus} -e '${cond}' ${SAM_NAME} | cut -f1  > ${NAMES}"""
	def FILTER_CMD = """LC_ALL=C fastq-namefilter ${NAMES} -"""

	"""
	set -euo pipefail
	${!meta.single_end ? ALIGNER_CMD_PE : ALIGNER_CMD_SE}
	mv Aligned.out.bam ${SAM_NAME}

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

def set_genome_size(n) {
	params.starGenomeSize = n
}
def ensure_genome_size(explicit_size, ref_genome_fa) {
	if (explicit_size && explicit_size > 0) {
		set_genome_size(explicit_size)
		return explicit_size // user provided a genome size
	} else if (ref_genome_fa) {
		set_genome_size(ref_genome_fa.size())
		return ref_genome_fa.size() // based on fasta file size
	} else {
		log.error "--genome_size must be provided\n${star_usage()}"
		exit 1
	}
}

def ensure_star_indexes(ref_indexes, ref_indexes_ercc,
						ref_genome, ref_genome_gtf, ercc, ercc_gtf) {
	if (ref_indexes && ref_indexes_ercc
		&& (indexes = file(ref_indexes)) && (indexes2 = file(ref_indexes_ercc))
		&& indexes.exists() && indexes2.exists()) {
	} else {
		(indexes, indexes2) = star_generate_indexes(file(ref_genome),
													file(ref_genome_gtf),
													file(ercc),
													file(ercc_gtf))
	}
	return [indexes, indexes2]
}

workflow {
	ensure_genome_size(params.genomeSize, file(params.refGenome))

	(indexes, indexes2) = ensure_star_indexes(params.refIndexes,
											  params.refIndexesErcc,
											  params.refGenome,
											  params.refGenomeGtf,
											  params.ercc,
											  params.erccGtf)

	fastqs = LOAD_METASHEET(params.metaIn)

	unmapped_reads = star(fastqs, indexes)

	SAVE_METASHEET(unmapped_reads, params.metaOut)
}
