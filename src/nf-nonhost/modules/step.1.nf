#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.accessions_list = "SRA_accession_list.test.txt"
params.parallel_downloads = 10
params.publish_dir = "$PWD"
params.publish_intermediate = true

params.sra_prefetch_options = ""
params.fastq_dump_options = ""

params.meta_out = 'step_1_sheet.csv'
include {
	SAVE_METASHEET;
} from './utils.nf'

process prefetch {
	debug true
	label 'sratools'
	maxForks params.parallel_downloads
	
	input:
	tuple val(meta), val(sra_id)
	
	output:
	tuple val(meta), path('[S,E,D]RR*[0-9]'), env(reads), env(sra_size)

	// fastq-dump wants sras in the current directory. this is a problem for
	// nf's usually directory-agnostic behavior - it could be that the input is
	// cached from a previous run and the absolute dir invisible,
	script:
	"""
	set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e
	prefetch --output-directory . --max-size 1t --force ALL \
		$params.sra_prefetch_options $sra_id
	vdb-validate -I no $sra_id 2> validate.txt
	vdb-dump --info $sra_id > info.txt
	sra_size=\$(awk -F': ' '/^size/{gsub(/,/,"",\$2);print \$2}' info.txt)
	reads=\$(awk -F': ' '/^SEQ/{gsub(/,/,"",\$2);print \$2}' info.txt)
	"""

	stub:
	"""
	mkdir -p ${sra_id}
	touch ${sra_id}/${sra_id}.sra
	"""
}

process fastq_dump {
	label 'sratools'

	publishDir "$params.publish_dir", enabled: params.publish_intermediate

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
		fasterq-dump --split-3 -e ${task.cpus} $params.fastq_dump_options --outdir fastq ${meta.id}
	
		# TODO: this doesn't work..
		rm -rf ${sra_file}
		"""
	else
		"""
		echo fasterq-dump encountered error, reverting to using fastq-dump
		set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e
		fastq-dump --split-3 --disable-multithreading --outdir fastq ${meta.id}

		rm -rf ${sra_file}
		"""

	stub:
	"""
	if [ -f ${meta.id}/${meta.id}.sra ]
	then touch fastq/${meta.id}.fastq; fi
	"""
}

// after the filter_barcodes step, meta must be annotated
// 'single_end: true' for single files.  future steps will
// use that meta value for conditional processing.
process filter_barcodes {
	label 'median'
	
	input:
	tuple val(meta), path(fastq)

	output:
	tuple val(meta), path('*.fastq.gz'), env(median), env(count), env(size)

	script:
	"""
	# one set of reads, or two?
	if [[ \$(ls -1 fastq/*.fastq | wc -l) -lt 2 ]]
	then # single
		IFS=\$'\\t' read -r -d \$'\\n' median differing count size <<< "\$(fastq-lengths summary fastq/${meta.id}.fastq)"
		echo $meta.id SE
		gzip -k6c fastq/${meta.id}.fastq > ${meta.id}.fastq.gz
	else # possibly paired, but may be scRNAseq barcodes
		IFS=\$'\\t' read -r -d \$'\\n' median1 differing1 count1 size1 <<< "\$(fastq-lengths summary fastq/${meta.id}_1.fastq)"
		IFS=\$'\\t' read -r -d \$'\\n' median differing count size <<< "\$(fastq-lengths summary fastq/${meta.id}_2.fastq)"
		if [ \$median1 -lt 32 ] && [ \$median -gt 80 ]
		then
			echo $meta.id scRNAseq
			rm fastq/${meta.id}_1.fastq # discard read 1 (cell barcode)
			mv fastq/${meta.id}_2.fastq fastq/${meta.id}.fastq
			gzip -k6c fastq/${meta.id}.fastq > ${meta.id}.fastq.gz
		else
			echo $meta.id PE
			gzip -k6c fastq/${meta.id}_1.fastq > ${meta.id}_1.fastq.gz & \
				gzip -k6c fastq/${meta.id}_2.fastq > ${meta.id}_2.fastq.gz
		fi
	fi
	gzip -t fastq/*.fastq.gz
	"""
}


workflow {
	accessions = channel.fromPath(params.accessions_list).splitText()
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

	SAVE_METASHEET(fastqs_2, params.meta_out)
}