#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.publishDir = "$PWD"
params.publishIntermediate = true

params.metaIn = 'step_1_sheet.csv'
params.metaOut = 'step_2_sheet.csv'
include {
	LOAD_METASHEET;
	SAVE_METASHEET;
} from './utils.nf'

process fastp {
	label 'fastp'

	input:
	tuple val(meta), path(fqs)

	output:
	tuple val(meta), path('*.trimmed.fastq'), emit: mates
	tuple val(meta), path("stats.txt"), emit: stats

	script:
	def extension = ".trimmed.fastq"
	def FASTP_CMD = """fastp --disable_quality_filtering --length_required 2 \
		--compression 6 --thread ${task.cpus} \
		--json fastp.json --html fastp.html """
	if (!meta.single_end) """
	trap 'echo "\$\$ Interrupt by external (OOM?), exiting."; exit 130' SIGINT

	gzip -dc ${fqs[0]} > fq1.fastq
	gzip -dc ${fqs[1]} > fq2.fastq
	${FASTP_CMD} \
		--detect_adapter_for_pe --in1 fq1.fastq --in2 fq2.fastq \
		--out1 ${meta.id}_1${extension}.staging \
		--out2 ${meta.id}_2${extension}.staging 2>stats.txt
	mv ${meta.id}_1${extension}.staging ${meta.id}_1${extension}
	mv ${meta.id}_2${extension}.staging ${meta.id}_2${extension}
	"""
	else if (meta.single_end) """
	trap 'echo "\$\$ Interrupt by external (OOM?), exiting."; exit 130' SIGINT

	gzip -dc ${fqs} > fq1.fastq
	${FASTP_CMD} \
		--in1 fq1.fastq --out1 ${meta.id}${extension}.staging 2>stats.txt
	mv ${meta.id}${extension}.staging ${meta.id}${extension}
	"""
}

process priceseqfilter {
	label 'price'

	publishDir "$params.publishDir/fastq/$meta.id", enabled: params.publishIntermediate
	
	input:
	tuple val(meta), path(fqs)

	output:
	tuple val(meta), path('*.PRICEfiltered.fastq.gz'), emit: mates
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
	
	${price_cmd} \
		-fp ${fqs} \
		-op ${meta.id}_1${extension} ${meta.id}_2${extension}
	gzip -nc ${meta.id}_1${extension} > ${meta.id}_1${extension}.gz.staging 
	gzip -nc ${meta.id}_2${extension} > ${meta.id}_2${extension}.gz.staging 
	mv ${meta.id}_1${extension}.gz.staging ${meta.id}_1${extension}.gz
	mv ${meta.id}_2${extension}.gz.staging ${meta.id}_2${extension}.gz
	gzip -t *${extension}.gz
	"""
	else if (meta.single_end)
	"""
	trap 'echo "\$\$ Interrupt by external (OOM?), exiting."; exit 130' SIGINT
	
	${price_cmd} \
		-f ${fqs} -o ${meta.id}${extension}
	gzip -nc ${meta.id}${extension} > ${meta.id}${extension}.gz.staging
	mv ${meta.id}${extension}.gz.staging ${meta.id}${extension}.gz
	gzip -t *${extension}.gz
	"""
}

workflow {
	fastqs = LOAD_METASHEET(params.metaIn)
	
	fastp(fastqs)
	priceseqfilter(fastp.out)
	
	SAVE_METASHEET(priceseqfilter.out, params.metaOut)
}
