#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.accessions_list = "SRA_accession_list.test.txt"
params.parallel_downloads = 10
params.publish_dir = "$PWD"
params.publish_intermediate = true
 
process prefetch {
	label 'sratools'
	label 'network_limited'
	maxForks params.parallel_downloads

	input:
	tuple val(meta), val(sra_id)
	
	output:
	tuple val(meta), path('[S,E,D]RR*[0-9]')
	
	// fastq-dump wants sras in the current directory. this is a problem for
	// nf's usually directory-agnostic behavior - it could be that the input is
	// cached from a previous run and the absolute dir invisible,
	script:
	"""
	set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e
	prefetch --output-directory . --max-size 1t --force ALL $sra_id
	"""
	// we can possibly parse total read count out of
	// vdb-dump --info /work/ERR1415401/ERR1415401.sra
	// SEQ    : 1,730,560

	stub:
	"""
	mkdir -p ${sra_id}
	touch ${sra_id}/${sra_id}.sra
	"""
}

process fastq_dump {
	label 'sratools'
	//	label 'mem_medium'
	cpus 8
	memory 30G
	publishDir "$params.publish_dir", enabled: params.publish_intermediate
	errorStrategy { task.exitStatus != 3 ? 'retry' : 'terminate' }
	maxRetries 1

	input:
	tuple val(meta), path(sra_file)

	output:
	tuple val(meta), path('fastq')

	beforeScript 'mkdir -p fastq'
	
	script:
	mem = task.memory.toString() - ~/ /
	if (task.attempt == 1)
		"""
		set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e
		fasterq-dump --split-3 -b 4G -c 4G -m 4G -e 8 --outdir fastq ${meta.id}
		gzip fastq/*.fastq

		# TODO: this doesn't work..
		rm -rf ${sra_file}
		"""
	else
		"""
		set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e
		echo fasterq-dump encountered error, reverting to using fastq-dump
		fastq-dump --split-3 --disable-multithreading --outdir fastq ${meta.id}
		gzip fastq/*.fastq
		rm -rf ${sra_file}
		"""

	stub:
	"""
	if [ -f ${meta.id}/${meta.id}.sra ]
	then touch fastq/${meta.id}.fastq.gz; fi
	"""
}

workflow {
	accessions = channel.fromPath(params.accessions_list).splitText()
		.map { acc -> [[id: acc.trim()], acc.trim()] }
	prefetch(accessions) | fastq_dump
}