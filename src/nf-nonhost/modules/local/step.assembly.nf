#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.fastqDir

params.publishDir = "$PWD"
params.publishIntermediate = true
params.genomeSize = null
params.maxMismatch = 0.3

params.metaIn = 'step_assembly_sheet.csv'
params.metaOut = 'step_assembly_sheet.csv'
include {
	LOAD_METASHEET;
	SAVE_METASHEET;
} from './utils.nf'


process spades {
	cpus 8
	memory "196GB"
	
	label 'spades'

	input:
	tuple val(meta), path(mategz)

	output:
	tuple val(meta), path(contigs)

	script:
	def SPADES_CMD = """spades.py --rna --threads ${task.cpus} \
		--memory ${task.memory.giga()} """
	if (!meta.single_end)
	"""
	${SPADES_CMD} -1 ${mategz[0]} -2 ${mategz[1]} -o ${meta.id}.staging
	mv ${meta.id}.staging ${meta.id}
	"""
}

workflow assemble {
	take: bioproject_mates

	bioproj_with_meta = bioproject_mates
		.map { meta, mates ->
			def new_meta = meta.clone()
			new_meta.single_end = mates.size() != 2
			[ new_meta, fastq ]
		}

	main:
	spades(bioproj_with_meta)
}

workflow {
	if (params.bioprojPath)
		assemble(
		Channel.fromPath("$params.bioprojPath/*", type: 'dir')
			.map { path ->
				def new_meta = [ id: path.getSimpleName(),
								fastq_size: files("$path/*.fastq")[0].size() ]
				[ new_meta, Channel.fromPath("$path/*.(fastq|fastq.gz)") ]
			})
}