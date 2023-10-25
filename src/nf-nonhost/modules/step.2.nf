#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.publish_dir = "$PWD"
params.publish_intermediate = true

params.fastp_options = ""
params.price_options = ""

params.meta_in = 'step_1_sheet.csv'
params.meta_out = 'step_2_sheet.csv'
include {
	LOAD_METASHEET;
	SAVE_METASHEET;
} from './utils.nf'

process fastp {
	label 'fastp'

	input:
	tuple val(meta), path(fqs)

	output:
	tuple val(meta), path('*.trimmed.fastq')

	script:
	def extension = ".trimmed.fastq"
	def FASTP_CMD = """fastp --disable_quality_filtering --length_required 2 \
		--compression 6 --thread ${task.cpus} \
		--json fastp.json --html fastp.html $params.fastp_options"""
	if (!meta.single_end)
	"""
	gzip -dc ${fqs[0]} > fq1.fastq
	gzip -dc ${fqs[1]} > fq2.fastq
	${FASTP_CMD} \
		--detect_adapter_for_pe --in1 fq1.fastq --in2 fq2.fastq \
		--out1 ${meta.id}_1${extension} --out2 ${meta.id}_2${extension} \
	rm ${fqs}
	"""
	else if (meta.single_end)
	"""
	gzip -dc ${fqs[0]} > fq1.fastq
	${FASTP_CMD} \
		--in1 fq1.fastq --out1 ${meta.id}${extension} \
	rm ${fqs}
	"""
}

process priceseqfilter {
	label 'price'

	publishDir "$params.publish_dir/fastq/$meta.id", enabled: params.publish_intermediate
	
	input:
	tuple val(meta), path(fqs)

	output:
	tuple val(meta), path('*.PRICEfiltered.fastq.gz')
	
// -a thread, -log c concise output, -fp input, -op output
// -rnf 90:    90 percentage of nucleotides in a read that must be called
// -rqf 85 0.98:       85% of the read length must be 98% accurate (parameterized after CZID, Joe uses 95% of the read length 98% accurate)
	script:
	def extension = ".PRICEfiltered.fastq"
	def price_cmd = """PriceSeqFilter -a ${task.cpus} -rnf 90 -rqf 85 0.98 \
		-log c $params.price_options"""
	if (!meta.single_end)
	"""
	${price_cmd} \
		-fp ${fqs} \
		-op ${meta.id}_1${extension} ${meta.id}_2${extension}
	gzip *${extension}
	rm ${fqs}
	"""
	else if (meta.single_end)
	"""
	${price_cmd} \
		-f ${fqs} -o ${meta.id}${extension}
	gzip *${extension}
	rm ${fqs}
	"""
}

workflow {
	fastqs = LOAD_METASHEET(params.meta_in)
	
	fastp(fastqs)
	priceseqfilter(fastp.out)
	
	SAVE_METASHEET(priceseqfilter.out, params.meta_out)
}
