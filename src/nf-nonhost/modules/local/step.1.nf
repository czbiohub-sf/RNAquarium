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


process prefetch {
	label 'sratools'
	maxForks params.parallelDownloads

	input:
	tuple val(meta), val(sra_id)

	output:
	tuple val(meta), path('[S,E,D]RR*[0-9]/*.sra*'), env(reads), env(sra_size), emit: sra
	tuple val(meta), path("info.txt"), emit: stats
	tuple val(meta), path("validate.txt"), emit: vdb_validate

	// fastq-dump wants sras in the current directory. this is a problem for
	// nf's usually directory-agnostic behavior - it could be that the input is
	// cached from a previous run and the absolute dir invisible,
	beforeScript = """module load mamba
	sleep \$((1 + RANDOM % 30))s"""

	script:
	"""
	trap 'echo "\$\$ Interrupt by external (OOM?), exiting."; exit 130' SIGINT

	set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e
	prefetch --output-directory staging --max-size 1t --force ALL $sra_id
	cd staging
	vdb-validate -I no $sra_id 2> ../validate.txt
	vdb-dump --info $sra_id > ../info.txt
	cd ..
	sra_size=\$(awk -F': ' '/^size/{gsub(/,/,"",\$2);print \$2}' info.txt)
	reads=\$(awk -F': ' '/^SEQ/{gsub(/,/,"",\$2);print \$2}' info.txt)
	trap -- '' SIGTERM
	mv staging/$sra_id $sra_id
	"""
}

process fastq_dump {
	label 'sratools'

	// omit 'fastq/$meta.id' because we output that entire folder here
	publishDir "$params.publishDir", enabled: params.publishIntermediate

	input:
	tuple val(meta), path("${meta.id}/${meta.id}.sra*")

	output:
	tuple val(meta), path("fastq/${meta.id}/*.fastq.gz"), emit: mates
	tuple val(meta), path("stats.txt"), emit: stats

	beforeScript = {"""module load mamba
		rm -rf \$NXF_SCRATCH || true
		local tmp_avail=\$(df -P "${params.tmp}" | tail -1 | awk '{print \$4}')
		tmp_choice="${params.backupTmp}"
		readarray -t ids < <(squeue -w \$SLURMD_NODENAME -o "%i" -S i | tail -n +2)
		PERJOBGUESS=67108864  # 64 GB? / 1024 (df reports 1K blocks)
		res_guess=0
		for (( i=0; i < \${#ids[@]}; i++ ));
		do
			if [[ \${ids[\$i]} -eq "\$SLURM_JOB_ID" ]];
			then
				res_guess=\$(( \$PERJOBGUESS * \$i ))
			fi
		done
		if [[ ${meta.sra_size} -ne 0 && \$(( ${meta.sra_size} * 50 * 2 )) -le \$(( (\$tmp_avail - \$res_guess) * 1024 )) ]]
		then
			tmp_choice="${params.tmp}"
		fi
		NXF_SCRATCH="\$(set +u; nxf_mktemp \$tmp_choice)"
					"""}


	script:
	mem = task.memory.toString() - ~/ /
	if (task.attempt == 1) """
	trap 'echo "\$\$ Interrupt by external (OOM?), exiting."; exit 130' SIGINT
	mkdir -p fastq
	set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e
	fasterq-dump --split-3 -e ${task.cpus} -m ${task.memory.toMega()-100}MB \
		--outdir fastq/${meta.id}.staging ${meta.id} 2>stats.txt

	trap -- '' SIGTERM
	${task.ext.gzipCmd} fastq/${meta.id}.staging/*.fastq
	mv fastq/${meta.id}.staging fastq/${meta.id}
	"""
	else """
	trap 'echo "\$\$ Interrupt by external (OOM?), exiting."; exit 130' SIGINT

	echo fasterq-dump encountered error, reverting to using fastq-dump
	mkdir -p fastq
	set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e
	fastq-dump --split-3 --disable-multithreading \
		--outdir fastq/${meta.id}.staging ${meta.id} 2>stats.txt

	trap -- '' SIGTERM
	${task.ext.gzipCmd} fastq/${meta.id}.staging/*.fastq
	mv fastq/${meta.id}.staging fastq/${meta.id}
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
	tuple val(meta), path("*.filtered.fastq.gz", arity: '1..2'), env(median), env(count), env(size)

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
	IFS=\$'\\t' read -r -d \$'\\n' median differing count size <<< "\$(fastq-lengths summary ${meta.id}.fastq)"
		echo $meta.id SE
		${task.ext.gzipCmd} -fnc ${meta.id}.fastq > ${meta.id}${SUFFIX}.staging
		mv ${meta.id}${SUFFIX}.staging ${meta.id}${SUFFIX}
	else # possibly paired, but may be scRNAseq barcodes
	IFS=\$'\\t' read -r -d \$'\\n' median1 differing1 count1 size1 <<< "\$(fastq-lengths summary ${meta.id}_1.fastq)"
	IFS=\$'\\t' read -r -d \$'\\n' median differing count size <<< "\$(fastq-lengths summary ${meta.id}_2.fastq)"
		if [ \$median1 -lt 32 ] && [ \$median -gt 80 ]
		then
			echo $meta.id scRNAseq
			rm ${meta.id}_1.fastq # discard read 1 (cell barcode)
			mv ${meta.id}_2.fastq ${meta.id}.fastq
			${task.ext.gzipCmd} -fnc ${meta.id}.fastq > ${meta.id}${SUFFIX}.staging
			mv ${meta.id}${SUFFIX}.staging ${meta.id}${SUFFIX}
		elif [ \$median1 -gt 80 ] && [ \$median -lt 32 ]
		then
			echo $meta.id scRNAseq
			rm ${meta.id}_2.fastq # discard read 2 (cell barcode)
			mv ${meta.id}_1.fastq ${meta.id}.fastq
			median=\$median1
			${task.ext.gzipCmd} -fnc ${meta.id}.fastq > ${meta.id}${SUFFIX}.staging
			mv ${meta.id}${SUFFIX}.staging ${meta.id}${SUFFIX}
		else
			echo $meta.id PE
			${task.ext.gzipCmd} -fnc ${meta.id}_1.fastq > ${meta.id}_1${SUFFIX}.staging
			${task.ext.gzipCmd} -fnc ${meta.id}_2.fastq > ${meta.id}_2${SUFFIX}.staging
			mv ${meta.id}_1${SUFFIX}.staging ${meta.id}_1${SUFFIX}
			mv ${meta.id}_2${SUFFIX}.staging ${meta.id}_2${SUFFIX}
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
