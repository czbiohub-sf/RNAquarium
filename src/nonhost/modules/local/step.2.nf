#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.publishDir = "$PWD"
params.publishIntermediate = true
params.cleanupScript = ""
params.tmp = null
params.backupTmp = null
params.backupScratchHack = false
params.nxfUnstageHack = false

params.extraAdapters = "extra-adapters.fasta"

params.metaIn = 'step_1_sheet.csv'
params.metaOut = 'step_2_sheet.csv'
include {
	LOAD_METASHEET;
	SAVE_METASHEET;
} from './utils.nf'

process fastp {
	label 'fastp'

	input:
	tuple val(meta), path(fqs, arity: '1..2')

	output:
	tuple val(meta), path('*.trimmed.fastq.gz', arity: '1..2'), env(fastp_after), emit: mates
	tuple val(meta), path("stats.txt"), emit: stats_txt
	tuple val(meta), path("fastp_stats.csv"), emit: stats_csv

	script:
	def extension = ".trimmed.fastq"
	def ex_adp = file(params.extraAdapters).exists() ? "--adapter_fasta ${params.extraAdapters}" : ""
	def FASTP_CMD = """fastp -q 17 -u 15 -n \$(( ${meta.readlen}/10 > 50 ? 50 : ${meta.readlen}/10 )) \
		--length_required 2 ${ex_adp} \
		--compression 6 --thread ${task.cpus} \
		--json fastp.json --html fastp.html """
	def STATS_CMD = """fastp_before=\$(sed -n '/Read1 before filtering:/{;n;p;}' stats.txt | cut -f3 -d' ')
	fastp_after=\$(sed -n '/Read1 after filtering:/{;n;p;}' stats.txt | cut -f3 -d' ')
	fastp_short=\$(sed -n '/reads failed due to too short:/{;p;}' stats.txt | cut -f7 -d' ')
	fastp_trimmed=\$(sed -n '/reads with adapter trimmed:/{;p;}' stats.txt | cut -f5 -d' ')
	printf "%s,%s,%s,%s," \$fastp_before \$fastp_after \$fastp_short \$fastp_trimmed > fastp_stats.csv """
	if (!meta.single_end) """
	trap 'echo "\$\$ Interrupt by external (OOM?), exiting."; exit 130' SIGINT

	${task.ext.gzipCmd} -fdc ${fqs[0]} > fq1.fastq &
	${task.ext.gzipCmd} -fdc ${fqs[1]} > fq2.fastq &
	wait
	${FASTP_CMD} \
		--detect_adapter_for_pe --in1 fq1.fastq --in2 fq2.fastq \
		--out1 ${meta.id}_1${extension}.staging \
		--out2 ${meta.id}_2${extension}.staging 2>stats.txt
	${task.ext.gzipCmd} -nc ${meta.id}_1${extension}.staging > ${meta.id}_1${extension}.gz.staging
	${task.ext.gzipCmd} -nc ${meta.id}_2${extension}.staging > ${meta.id}_2${extension}.gz.staging
	mv ${meta.id}_1${extension}.gz.staging ${meta.id}_1${extension}.gz
	mv ${meta.id}_2${extension}.gz.staging ${meta.id}_2${extension}.gz

	${STATS_CMD}
	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
	else if (meta.single_end) """
	trap 'echo "\$\$ Interrupt by external (OOM?), exiting."; exit 130' SIGINT

	${task.ext.gzipCmd} -fdc ${fqs} > fq1.fastq
	${FASTP_CMD} \
		--in1 fq1.fastq --out1 ${meta.id}${extension}.staging 2>stats.txt
	${task.ext.gzipCmd} -nc ${meta.id}${extension}.staging > ${meta.id}${extension}.gz.staging
	mv ${meta.id}${extension}.gz.staging ${meta.id}${extension}.gz

	${STATS_CMD}
	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
}

process priceseqfilter {
	label 'price'

	publishDir "$params.publishDir/fastq/$meta.id", enabled: params.publishIntermediate

	input:
	tuple val(meta), path(fqs, arity: '1..2')

	output:
	tuple val(meta), path('*.PRICEfiltered.fastq.gz', arity: '1..2'), emit: mates
	tuple val(meta), stdout, emit: stats

// -a thread, -log c concise output, -fp input, -op output
// -rnf 90:    90 percentage of nucleotides in a read that must be called
// -rqf 85 0.98:       85% of the read length must be 98% accurate (parameterized after CZID, Joe uses 95% of the read length 98% accurate)
	script:
	def extension = ".PRICEfiltered.fastq"
	def price_cmd = """PriceSeqFilter -a ${task.cpus} -rnf 90 -rqf 85 0.98 \
		-log c """
	if (!meta.single_end)
	"""
	trap 'echo "\$\$ Interrupt by external (OOM?), exiting."; exit 130' SIGINT
	gunzip -dc ${fqs[0]} > fq1.fastq
	gunzip -dc ${fqs[1]} > fq2.fastq

	${price_cmd} \
		-fp fq1.fastq fq2.fastq \
		-op ${meta.id}_1${extension} ${meta.id}_2${extension} | tee stats.txt
	${task.ext.gzipCmd} -nc ${meta.id}_1${extension} > ${meta.id}_1${extension}.gz.staging
	${task.ext.gzipCmd} -nc ${meta.id}_2${extension} > ${meta.id}_2${extension}.gz.staging
	mv ${meta.id}_1${extension}.gz.staging ${meta.id}_1${extension}.gz
	mv ${meta.id}_2${extension}.gz.staging ${meta.id}_2${extension}.gz
	${task.ext.gzipCmd} -t *${extension}.gz

	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
	else if (meta.single_end)
	"""
	trap 'echo "\$\$ Interrupt by external (OOM?), exiting."; exit 130' SIGINT
	gunzip -dc ${fqs} > fq1.fastq

	${price_cmd} \
		-f fq1.fastq -o ${meta.id}${extension} | tee stats.txt
	${task.ext.gzipCmd} -nc ${meta.id}${extension} > ${meta.id}${extension}.gz.staging
	mv ${meta.id}${extension}.gz.staging ${meta.id}${extension}.gz
	${task.ext.gzipCmd} -t *${extension}.gz

	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
}

workflow {
	fastqs = LOAD_METASHEET(params.metaIn)

	fastp(fastqs)
	priceseqfilter(fastp.out)

	SAVE_METASHEET(priceseqfilter.out, params.metaOut)
}
