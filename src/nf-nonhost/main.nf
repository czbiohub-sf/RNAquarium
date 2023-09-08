#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.help = false
params.h = false
params.accessions_list = ""
params.parallel_downloads = 10
params.publish_dir = "$PWD"
params.publish_intermediate = false

include {
	prefetch;
	fastq_dump;
} from './step.1.nf' params(
	parallel_downloads: params.parallel_downloads,
	publish_dir: params.publish_dir,
	publish_intermediate: params.publish_intermediate
)
include {
	filter_barcodes;
	fastp;
	priceseqfilter;
} from './step.2.nf' params(
	publish_dir: params.publish_dir,
	publish_intermediate: params.publish_intermediate
)

//dependencies:
// nextflow seems to have two easy options
// - processes can define containers to run in
// - programs in the ./bin/ directory are available (?) to processes
// so we may want to allow running in either 'use containers' or
// 'i have pre-installed the dependencies' mode.
// this gets a little awkward with our hpc, where the preferred method
// of using some software is `module load`.
// so for specific cases it may be a 
// 'i have pre-installed the dependencies without the provided instructions
// and don't hold the pipeline accountable for issues with that configuration'
// right, perhaps 'instructions to install to the bin/ path that nextflow wants
// htseq is annoying!!!!!!


def helpMessage() {
	log.info """
	--accessions_list path    file containing sra accessions to process, one per line
	                            (required)
	--parallel_prefetch n     maximum sra prefetch downloads to run at once
	                            (default: 10)
	--publish_dir path        path to write useful intermediate output of each step to
	--help, -h                print this text and exit
	
	-with-docker              use docker containers to run commands
	-with-singularity         use singularity containers to run commands
                                (avoid most pre-installation procedure)
	"""
}

workflow {
	if (params.help || params.h || !params.accessions_list) {
		helpMessage()
		exit params.help || params.h ? 0 : 1
	}
	
	main:
	// step 1
	accessions = channel.fromPath(params.accessions_list).splitText()
	fastqs = accessions | prefetch | fastq_dump

	// step 2
	filter_barcodes(fastqs) | fastp | priceseqfilter

}