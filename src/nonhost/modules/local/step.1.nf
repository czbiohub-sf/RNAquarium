#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.accessionsList = "SRA_accession_list.test.txt"
params.parallelDownloads = 10
params.publishDir = "$PWD"
params.publishIntermediate = true
params.cleanupScript = ""
params.tmp = null
params.backupTmp = null
params.backupScratchHack = false
params.nxfUnstageHack = false

params.metaOut = 'step_1_sheet.csv'
include {
	SAVE_METASHEET;
} from './utils.nf'


process download {
	label 'sratools'
	maxForks params.parallelDownloads

	input:
	tuple val(meta), val(sra_id)

	output:
	tuple val(meta), path("fastq/${sra_id}/*.fastq.gz", arity: '1..2'), env(reads), env(sra_size), env(median1), env(median2), env(count1), env(size1), emit: mates
	tuple val(meta), path("info.txt"), emit: srastats
	tuple val(meta), path("stats.txt"), emit: dumpstats

	beforeScript = {"""${task.ext.extraBeforeScript ?: ""}
					sleep \$((1 + RANDOM % 30))s"""}

	script:
	def PROLOGUE = """
	trap 'echo "\$\$ Interrupt by external (OOM?), exiting."; exit 130' SIGINT
	set -v
	
	# prefetch won't run without config, we can't control much, but at least initialize it
	set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e
	mkdir -p fastq
	"""
	def PREFETCH = """
	# we get sralite most of the time, but don't want to fail if it doesn't exist
	prefetch --output-directory staging --max-size 1t --force ALL ${sra_id}
	cd staging
	vdb-dump --info ${sra_id} > ../info.txt
	cd ..
	sra_size=\$(awk -F': ' '/^size/{gsub(/,/,"",\$2);print \$2}' info.txt)
	reads=\$(awk -F': ' '/^SEQ/{gsub(/,/,"",\$2);print \$2}' info.txt)
	"""
	// consider collecting vdb-dump $sra_id -R 1-100000 -C READ_TYPE,READ_LEN -ftab | \
	// sed "s/\(SRA_READ_TYPE_\|IOLOGICAL\|ECHNICAL\|FORWARD\|REVERSE\||\)//g" | \
	// sort --parallel=4 | uniq -c | sort -rh
	// as well
	def EPILOGUE = """
	trap -- '' SIGTERM
	"""
	def FASTQ_DUMP = (task.attempt <= 2) ? """
	fasterq-dump staging/${sra_id} --split-3 --temp /dev/shm -x -e ${task.cpus} \
		-b 4M -c 32M \
		-m ${4096*task.attempt}M \
		--seq-defline '@\$ac.\$si/\$ri' --qual-defline '+' \
		--outdir fastq/${sra_id}.staging 2>stats.txt
	""" : """
	echo fasterq-dump encountered error, reverting to using fastq-dump
	fastq-dump staging/${sra_id} --split-3 \
		--defline-seq '@\$ac.\$si/\$ri' --defline-qual '+' \
		--outdir fastq/${sra_id}.staging 2>stats.txt
	"""
	"""
	${PROLOGUE}
	${PREFETCH}
	${FASTQ_DUMP}
	rm -rf ./staging/${sra_id}

	# filter obvious barcode reads
	f=fastq/${sra_id}.staging/${sra_id}.fastq
	f1=fastq/${sra_id}.staging/${sra_id}_1.fastq
	f2=fastq/${sra_id}.staging/${sra_id}_2.fastq
	if [[ ! ( -e \${f1} && -e \${f2} ) ]]; then  # single read file
		IFS=\$'\\t' read -rd \$'\\n' median1 differing1 count1 size1 <<<"\$(fastq-lengths summary \${f})"
		echo $meta.id SE
		rm -f \${f1} \${f2}
	else # possibly paired, but may be scRNAseq barcodes
		IFS=\$'\\t' read -rd \$'\\n' median1 differing1 count1 size1 <<< "\$(fastq-lengths summary 2000000 \${f1})"
		IFS=\$'\\t' read -rd \$'\\n' median2 differing2 count2 size2 <<< "\$(fastq-lengths summary 2000000 \${f2})"
		if [ \$median1 -lt 32 ] && [ \$median2 -gt 80 ]; then
			echo $meta.id scRNAseq
			rm \${f1} # discard read 1 (cell barcode)
			rm -f \${f}
			mv \${f2} \${f}
			size1=\$size2
		elif [ \$median1 -gt 80 ] && [ \$median2 -lt 32 ]; then
			echo $meta.id scRNAseq
			rm \${f2} # discard read 2 (cell barcode)
			rm -f \${f}
			mv \${f1} \${f}
		else
			echo $meta.id PE
			rm -f \${f}
			size1=\$(bc<<<"\$size2+\$size1")
		fi
	fi
	
	${task.ext.gzipCmd} fastq/${sra_id}.staging/*.fastq
	mv fastq/${sra_id}.staging fastq/${sra_id}
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
	// we assume that $fastqs are named {meta.id}_1.fastq, {meta.id}_2.fastq, {meta.id}.fastq only.
	tuple val(meta), path(fastqs)

	output:
	tuple val(meta), path("*.filtered.fastq.gz", arity: '1..2'), env(median1), env(median2), env(count1), env(size1)

	script:
	def SUFFIX = ".filtered.fastq.gz"
	"""
	trap 'echo "\$\$ Interrupt by external (OOM?), exiting."; exit 130' SIGINT

	for file in $fastqs; do
		${task.ext.gzipCmd} -dc \$file > \${file%.*}
	done

	# one set of reads, or two?
	if [[ ! ( -e ${meta.id}_1.fastq && -e ${meta.id}_2.fastq ) ]]
	then # single
		IFS=\$'\\t' read -r -d \$'\\n' median1 differing1 count1 size1 <<< "\$(fastq-lengths summary ${meta.id}.fastq)"
		echo $meta.id SE
		${task.ext.gzipCmd} -fnc ${meta.id}.fastq > ${meta.id}${SUFFIX}.staging
		mv ${meta.id}${SUFFIX}.staging ${meta.id}${SUFFIX}
	else # possibly paired, but may be scRNAseq barcodes
		IFS=\$'\\t' read -r -d \$'\\n' median1 differing1 count1 size1 <<< "\$(fastq-lengths summary ${meta.id}_1.fastq)"
		IFS=\$'\\t' read -r -d \$'\\n' median2 differing2 count2 size2 <<< "\$(fastq-lengths summary ${meta.id}_2.fastq)"
		if [ \$median1 -lt 32 ] && [ \$median2 -gt 80 ]; then
			echo $meta.id scRNAseq
			rm ${meta.id}_1.fastq # discard read 1 (cell barcode)
			mv ${meta.id}_2.fastq ${meta.id}.fastq
			${task.ext.gzipCmd} -fnc ${meta.id}.fastq > ${meta.id}${SUFFIX}.staging
			mv ${meta.id}${SUFFIX}.staging ${meta.id}${SUFFIX}
			size1=\$size2
		elif [ \$median1 -gt 80 ] && [ \$median2 -lt 32 ]; then
			echo $meta.id scRNAseq
			rm ${meta.id}_2.fastq # discard read 2 (cell barcode)
			mv ${meta.id}_1.fastq ${meta.id}.fastq
			${task.ext.gzipCmd} -fnc ${meta.id}.fastq > ${meta.id}${SUFFIX}.staging
			mv ${meta.id}${SUFFIX}.staging ${meta.id}${SUFFIX}
		else
			echo $meta.id PE
			${task.ext.gzipCmd} -fnc ${meta.id}_1.fastq > ${meta.id}_1${SUFFIX}.staging
			${task.ext.gzipCmd} -fnc ${meta.id}_2.fastq > ${meta.id}_2${SUFFIX}.staging
			mv ${meta.id}_1${SUFFIX}.staging ${meta.id}_1${SUFFIX}
			mv ${meta.id}_2${SUFFIX}.staging ${meta.id}_2${SUFFIX}
			size1=\$(bc<<<"\$size2+\$size1")
		fi
	fi
	${task.ext.gzipCmd} -t *${SUFFIX}

	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
}


workflow {
	accessions = channel.fromPath(params.accessionsList).splitText()
		.map { acc -> [[id: acc.trim()], acc.trim()] }

	fastqs_2 = download(accessions).mates
		.map { meta, fastq, reads, sra_size, median1, median2, count, fsize ->
			def new_meta = [id: meta.id,
							reads: count.toLong(),
							sra_size: sra_size.toLong(),
							readlen: median1.toLong(),
							readlen_2: median2 != "" ? median2.toLong() : "",
							fastq_size: fsize.toLong(),
							single_end: fastq.size() != 2,
							cleanup: "",
							cleanup_later: "${sra.toString()}"]
			[ new_meta, sra ]
		}


	SAVE_METASHEET(fastqs_2, params.metaOut)
}
