#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.bioprojPath = null
params.publishDir = "$PWD"
params.publishIntermediate = true
params.genomeSize = null
params.tmp = null
params.backupTmp = null
params.backupScratchHack = false
params.nxfUnstageHack = false

params.metaIn = 'step_assembly_sheet.csv'
params.metaOut = 'step_assembly_sheet.csv'
include {
	LOAD_METASHEET;
	SAVE_METASHEET;
} from './utils.nf'


process spades {
    label 'spades'

	debug true
	publishDir "$params.publishDir/SPAdes/"
	cpus 8
	memory "196GB"
	
	label 'spades'

	input:
	tuple val(meta), path(mategz)

	output:
    tuple val(meta), path("$meta.id")
	
	script:
    def SPADES_CMD = """spades.py --rnaviral --threads ${task.cpus} \
		--memory ${task.memory.toGiga()} """
	if (!meta.single_end)
	"""
	${SPADES_CMD} -1 ${mategz[0]} -2 ${mategz[1]} -o ${meta.id}.staging
	mv ${meta.id}.staging ${meta.id}
	"""
}

workflow assemble {
    take: bioprojects_dir

	main:
	log.info bioprojects_dir
	bioproject_mates = Channel.fromPath(file(bioprojects_dir).resolve('*'), followLinks: true)
	    .view()
	    .map { path ->
			def new_meta = [ id: path.getSimpleName(),
							fastq_size: files("$path/*.fastq.gz")[0].size() ]
			[ new_meta, Channel.fromPath("$path/*.fastq.gz") ]
		}
	bioproject_mates.view()
	
	bioproj_with_meta = bioproject_mates
		.map { meta, mates ->
			def new_meta = meta.clone()
			new_meta.single_end = mates.size() != 2
			[ new_meta, mates ]
		}

	spades(bioproj_with_meta)
}

workflow {
    if (params.bioprojPath) {
	    bioprojects_dir = file(params.bioprojPath)
	    bioproject_mates = Channel.fromPath("$bioprojects_dir/*", type: 'dir', followLinks: true)
	        .map { path ->
			    def new_meta = [ id: path.getSimpleName(),
				                fastq_size: files("$path/*.fastq.gz")[0].size() ]
			    [ new_meta, files("$path/*.fastq.gz") ]
			}
	    bioproj_with_meta = bioproject_mates
	        .map { meta, mates ->
			    def new_meta = meta.clone()
			    new_meta.single_end = mates.size() != 2
			    [ new_meta, mates ]
			}

	    spades(bioproj_with_meta)
	}
}
