#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.accessions_list = "SRA_accession_list.test.txt"

process prefetch {
	debug true
//	cache true
//	cpus 4
//	memory 24G
	// todo: scratch directive for hpc
	// todo: container 'fastq_dump'
	input:
	val sra_id

	output:
	path '[S,E,D]RR*[0-9]'

	//beforeScript 'mkdir -p prefetched'
	
	// fastq-dump wants sras in the current directory. this is a problem for
	// nf's usually directory-agnostic behavior - it could be that the input is
	// cached from a previous run and the absolute dir invisible,
	script:
	out_file = sra_id.trim()
	"""
	prefetch --output-directory . --max-size 1t --force ALL $sra_id 
	"""
}

process fastq_dump {
	debug true
	//memory '24GB'
	memory '16GB'

	input:
	path sra_id
	
	output:
	path 'fastq'

	beforeScript 'mkdir -p fastq'
	
	script:
	if (task.attempt == 1)
		"""
		fasterq-dump --split-3 --mem ${task.memory.toString().replace(" ", "")} --outdir fastq ${sra_id.baseName}
		gzip fastq/*.fastq
		"""
	else
		"""
		echo fasterq-dump encountered error, reverting to using fastq-dump
		fastq-dump --split-3 --disable-multithreading --outdir fastq ${sra_id.baseName}
		gzip fastq/*.fastq
		"""

	stub:
	"""
	touch fastq/${sra_file.baseName}.fastq.gz
	"""
}

workflow {
	accessions = channel.fromPath(params.accessions_list).splitText()
	prefetch(accessions) | fastq_dump
}