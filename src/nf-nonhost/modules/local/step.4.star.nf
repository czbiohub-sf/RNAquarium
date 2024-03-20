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
	publishDir "$params.publishDir/STAR_out/$meta.id", enabled: params.publishIntermediate
	
	input:
	tuple val(meta), path(fqgz, arity: '1..2')
	path indexes_dir // "Danio_rerio.GRCz11.108.ERCC"

	output:
	tuple val(meta), path("?E/Unmapped.out.mate?.gz", arity: '1..2'), emit: mates
	tuple val(meta), path("Log.final.out"), emit: stats

	script:
	def loadType = params.starUseSharedMem ? "LoadAndRemove" : "NoSharedMemory"
	def STAR_CMD = """STAR --outFilterMultimapNmax 99999 --outFilterMismatchNmax 999 \
		--outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5 \
		--outSAMmode None \
		--clip3pNbases 0 --limitOutSJcollapsed 200000000 \
		--genomeLoad $loadType \
		--outReadsUnmapped Fastx \
		--runThreadN ${task.cpus} """
	if (!meta.single_end)
	"""
	${STAR_CMD} \
		--genomeDir ${indexes_dir} \
		--readFilesCommand ${task.ext.gzipCmd} -dc \
		--outFileNamePrefix PE/ \
		--readFilesIn ${fqgz[0]} ${fqgz[1]}

	${task.ext.gzipCmd} -n PE/Unmapped.out.mate1
	${task.ext.gzipCmd} -n PE/Unmapped.out.mate2
	${task.ext.gzipCmd} -t PE/Unmapped.out.mate*.gz
	mv PE/Log.final.out Log.final.out

	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
	else if (meta.single_end)
	"""
	${STAR_CMD} \
		--genomeDir ${indexes_dir} \
		--readFilesCommand ${task.ext.gzipCmd} -dc \
		--outFileNamePrefix SE/ \
		--readFilesIn ${fqgz}

	${task.ext.gzipCmd} -n SE/Unmapped.out.mate1
	${task.ext.gzipCmd} -t SE/Unmapped.out.mate1.gz
	mv SE/Log.final.out Log.final.out

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
