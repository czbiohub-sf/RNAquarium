#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def star_usage() {
	return """STAR index generation parameters (one of the below options is required):
Option 1:
--ref_indexes_ercc path   pre-generated STAR indexes, with ERCC spike-in controls
--ref_indexes path        pre-generated STAR indexes for reference genome, without spike-in
--genome_size n           reference genome size, in bytes. (fa file size is good enough)
Option 2:
--ref_genome path         path to ensembl reference genome fasta
	e.g. Danio_rerio.GRCz11.dna_sm.primary_assembly.fa
--ref_genome_gtf path     path to ensembl reference genome annotations
    e.g. Danio_rerio.GRCz11.108.gtf
--ercc_fa path            path to ERCC spike-in control sequences
                            (default: ERCC92.fa)
--ercc_gtf path           path to ERCC spike-in control annotations
                            (default: ERCC92.gtf)
"""
}

params.publishDir = "$PWD"
params.publishIntermediate = true
params.genomeSize = null
params.refIndexesErcc = null // "Danio_rerio.GRCz11.108.ERCC"
params.refIndexes = null //"Danio_rerio.GRCz11.108"

params.refGenome = "Danio_rerio.GRCz11.dna_sm.primary_assembly.fa"
params.refGenomeGtf = "Danio_rerio.GRCz11.108.gtf"
params.erccFa = "ERCC92.fa"
params.erccGtf = "ERCC92.gtf"

params.starOptions = ""
params.starIndexGenOptions = ""

include {
	star_generate_indexes;
} from './step.0.generate_indexes.nf' params(
	publishDir: params.publishDir,
	starIndexGenOptions: params.starIndexGenOptions
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
	tuple val(meta), path(fqgz)
	path indexes_dir // "Danio_rerio.GRCz11.108.ERCC"

	output:
	tuple val(meta), path("?E/Unmapped.out.mate?.gz")

	script:
	def STAR_CMD = """STAR --outFilterMultimapNmax 99999 --outFilterMismatchNmax 999 \
		--outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5 \
		--outSAMmode None \
		--clip3pNbases 0 --limitOutSJcollapsed 200000000 \
		--genomeLoad NoSharedMemory --outReadsUnmapped Fastx \
		--runThreadN ${task.cpus} \
		${params.starOptions}"""
	if (!meta.single_end)
	"""
	${STAR_CMD} \
		--genomeDir ${indexes_dir} \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix PE/ \
		--readFilesIn ${fqgz[0]} ${fqgz[1]}

	gzip PE/Unmapped.out.mate*
	gzip -t PE/Unmapped.out.mate*.gz
	"""
	else if (meta.single_end)
	"""
	${STAR_CMD} \
		--genomeDir ${indexes_dir} \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix SE/ \
		--readFilesIn ${fqgz}

	gzip SE/Unmapped.out.mate1
	gzip -t SE/Unmapped.out.mate1.gz
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
	if (ref_indexes && ref_indexes_ercc) {
		(indexes, indexes2) = [file(ref_indexes_ercc), file(ref_indexes)]
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
