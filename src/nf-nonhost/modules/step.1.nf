#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.accessionsList = "SRA_accession_list.test.txt"
params.parallelDownloads = 10
params.publishDir = "$PWD"
params.publishIntermediate = true

params.sraPrefetchOptions = ""
params.fastqDumpOptions = ""

params.metaOut = 'step_1_sheet.csv'
include {
	SAVE_METASHEET;
} from './utils.nf'

process prefetch {
	debug true
	label 'sratools'
	maxForks params.parallelDownloads

	input:
	tuple val(meta), val(sra_id)

	output:
	tuple val(meta), path('[S,E,D]RR*[0-9]'), env(reads), env(sra_size), emit: sra
	tuple val(meta), path("info.txt"), emit: stats

	// fastq-dump wants sras in the current directory. this is a problem for
	// nf's usually directory-agnostic behavior - it could be that the input is
	// cached from a previous run and the absolute dir invisible,
	script:
	"""
	set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e
	prefetch --output-directory staging --max-size 1t --force ALL \
		$params.sraPrefetchOptions $sra_id
	cd staging
	vdb-validate -I no $sra_id 2> ../validate.txt
	vdb-dump --info $sra_id > ../info.txt
	cd ..
	sra_size=\$(awk -F': ' '/^size/{gsub(/,/,"",\$2);print \$2}' info.txt)
	reads=\$(awk -F': ' '/^SEQ/{gsub(/,/,"",\$2);print \$2}' info.txt)
	trap -- '' SIGTERM
	mv staging/$sra_id $sra_id
	"""

	stub:
	"""
	mkdir -p ${sra_id}
	touch ${sra_id}/${sra_id}.sra
	"""
}

process fastq_dump {
	label 'sratools'

	// omit 'fastq/$meta.id' because we output that entire folder here
	publishDir "$params.publishDir", enabled: params.publishIntermediate

	input:
	tuple val(meta), path(sra_file)

	output:
	tuple val(meta), path("fastq/${meta.id}/*.fastq"), emit: mates
	tuple val(meta), path("stats.txt"), emit: stats

	beforeScript 'mkdir -p fastq'

	script:
	mem = task.memory.toString() - ~/ /
	if (task.attempt == 1) """
	set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e
	fasterq-dump --split-3 -e ${task.cpus} $params.fastqDumpOptions \
		--outdir fastq/${meta.id}.staging ${meta.id} 2>stats.txt
	
	# TODO: this doesn't work..
	rm -rf ${sra_file}

	trap -- '' SIGTERM
	mv fastq/${meta.id}.staging fastq/${meta.id}
	"""
	else """
	echo fasterq-dump encountered error, reverting to using fastq-dump
	set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e
	fastq-dump --split-3 --disable-multithreading --outdir fastq/${meta.id}.staging ${meta.id}
	rm -rf ${sra_file}

	trap -- '' SIGTERM
	mv fastq/${meta.id}.staging fastq/${meta.id}
	"""

	stub:
	"""
	if [ -f ${meta.id}/${meta.id}.sra ]
	then touch fastq/${meta.id}.fastq; fi
	"""
}

process check_direct_fastqs {
	input:
	tuple val(meta), path(fastqs)

	output:
	tuple val(meta), path("$fastqs/*.fastq")

	script:
	"""
	"""
}

// after the filter_barcodes step, meta must be annotated
// 'single_end: true' for single files.  future steps will
// use that meta value for conditional processing.
process filter_barcodes {
	label 'median'

	input:
	tuple val(meta), path(fastqs)

	output:
	tuple val(meta), path('*.fastq.gz'), env(median), env(count), env(size)

	script: """
	# one set of reads, or two?
	if [ \$(echo "$fastqs" | wc -w) -ne 2 ]
	then # single
	IFS=\$'\\t' read -r -d \$'\\n' median differing count size <<< "\$(fastq-lengths summary $fastqs)"
		echo $meta.id SE
		gzip -k6c $fastqs > ${meta.id}.fastq.gz.staging
		mv ${meta.id}.fastq.gz.staging ${meta.id}.fastq.gz
	else # possibly paired, but may be scRNAseq barcodes
	IFS=\$'\\t' read -r -d \$'\\n' median1 differing1 count1 size1 <<< "\$(fastq-lengths summary ${fastqs[0]})"
	IFS=\$'\\t' read -r -d \$'\\n' median differing count size <<< "\$(fastq-lengths summary ${fastqs[1]})"
		if [ \$median1 -lt 32 ] && [ \$median -gt 80 ]
		then
			echo $meta.id scRNAseq
			rm ${fastqs[0]} # discard read 1 (cell barcode)
			mv ${fastqs[1]} ${meta.id}.fastq
			gzip -k6c ${meta.id}.fastq > ${meta.id}.fastq.gz.staging
			mv ${meta.id}.fastq.gz.staging ${meta.id}.fastq.gz
		else
			echo $meta.id PE
			gzip -k6c ${fastqs[0]} > ${meta.id}_1.fastq.gz.staging
			gzip -k6c ${fastqs[1]} > ${meta.id}_2.fastq.gz.staging
			mv ${meta.id}_1.fastq.gz.staging ${meta.id}_1.fastq.gz
			mv ${meta.id}_2.fastq.gz.staging ${meta.id}_2.fastq.gz
		fi
	fi
	gzip -t *.fastq.gz
	"""
}


workflow {
	accessions = channel.fromPath(params.accessionsList).splitText()
		.map { acc -> [[id: acc.trim()], acc.trim()] }

	sra = prefetch(accessions)
		.map { meta, sra, reads, sra_size ->
			def new_meta = [id: meta.id,
							reads: reads.toInteger(),
							sra_size: sra_size.toInteger() ]
			[ new_meta, sra ]
		}

	fastqs = fastq_dump(sra)

	fastqs_2 = filter_barcodes(fastqs)
		.map { meta, fastq, median, count, fsize ->
			def new_meta = meta.clone()
			new_meta.reads = count.toInteger()
			new_meta.readlen = median.toInteger()
			new_meta.fastq_size = fsize.toInteger()
			new_meta.single_end = fastq.size() != 2
			[ new_meta, fastq ]
		}

	SAVE_METASHEET(fastqs_2, params.metaOut)
}