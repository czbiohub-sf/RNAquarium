#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.accessions_list = "SRA_accession_list.test.txt"
params.parallel_downloads = 10
params.publish_dir = "$PWD"
params.publish_intermediate = true
/*
#SBATCH --job-name=fastqDump
#SBATCH --time=14-00:00:00
#SBATCH --array=1-10%10
#SBATCH --partition cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH -e slurm.out/slurm-%A_%a.err
#SBATCH -o slurm.out/slurm-%A_%a.out
*/
 
process prefetch {
	label 'sratools'
	maxForks params.parallel_downloads

	input:
	val sra_id

	output:
	path '[S,E,D]RR*[0-9]'
	
	// fastq-dump wants sras in the current directory. this is a problem for
	// nf's usually directory-agnostic behavior - it could be that the input is
	// cached from a previous run and the absolute dir invisible,
	script:
	"""
	set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e
	prefetch --output-directory . --max-size 1t --force ALL $sra_id 
	"""

	stub:
	"""
	touch $sra_id/$sra_id.sra
	"""
}

process fastq_dump {
	label 'sratools'
	publishDir params.publish_dir, enabled: params.publish_intermediate
	errorStrategy { task.exitStatus != 3 ? 'retry' : 'terminate' }
	maxRetries 1

	memory '16GB'

	input:
	path sra_file
	
	output:
	path 'fastq'

	beforeScript 'mkdir -p fastq'
	
	script:
	mem = task.memory.toString() - ~/ /
	if (task.attempt == 1)
		"""
		set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e
		fasterq-dump --split-3 --mem $mem --outdir fastq ${sra_file.baseName}
		gzip fastq/*.fastq

		# TODO: this doesn't work..
		rm -rf ${sra_file.baseName}/
		"""
	else
		"""
		set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e
		echo fasterq-dump encountered error, reverting to using fastq-dump
		fastq-dump --split-3 --disable-multithreading --outdir fastq ${sra_file.baseName}
		gzip fastq/*.fastq
		rm -rf ${sra_file.baseName}/
		"""

	stub:
	"""
	if [ -f ${sra_file.baseName}/${sra_file.baseName}.sra ]
	then touch fastq/${sra_file.baseName}.fastq.gz; fi
	"""
}

workflow {
	accessions = channel.fromPath(params.accessions_list).splitText()
	prefetch(accessions) | fastq_dump
}