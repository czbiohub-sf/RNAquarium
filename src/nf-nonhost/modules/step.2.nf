#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.fastq_path = "fastq/*/"
params.publish_dir = "$PWD"
params.publish_intermediate = true

// after the filter_barcodes step, meta must be annotated
// 'single_end: true' for single files.  future steps will
// use that meta value for conditional processing.
process filter_barcodes {
	label 'median'
	label 'mem_low'

	input:
	tuple val(meta), path(fastq)

	output:
	tuple val(meta), path('fastq/*.fastq')

	script:
	"""
	# unzipped may be cached
	if [ \$(ls fastq/*.gz | wc -w) -gt 0 ]; then gunzip -q fastq/*.gz; fi

	# one set of reads, or two?
	if [[ \$(ls -1 fastq | wc -l) -lt 2 ]]
	then # single
		#echo \$(fastq-lengths median 1 fastq/${meta.id}.fastq)
		echo $meta.id SE
	else # possibly paired, but may be scRNAseq barcodes
		if [ \$(fastq-lengths median fastq/${meta.id}_1.fastq) -lt 32 ] && \
			[ \$(fastq-lengths median fastq/${meta.id}_2.fastq) -gt 80 ]
		then
			echo $meta.id scRNAseq
			rm fastq/${meta.id}_1.fastq # discard read 1 (cell barcode)
			mv fastq/${meta.id}_2.fastq fastq/${meta.id}.fastq
		else
			echo $meta.id PE
		fi
	fi
	"""
}

process fastp {
	label 'fastp'
	//label 'mem_medium'
	cpus 8
	memory 8.GB

	input:
	tuple val(meta), path(fqs)

	output:
	tuple val(meta), path('*.trimmed.fastq')

	script:
	def extension = ".trimmed.fastq"
	def fastp_cmd = """fastp --disable_quality_filtering --disable_length_filtering \
		--compression 6 --thread 8 \
		--json fastp.json --html fastp.html"""
	if (!meta.single_end)
	"""
	${fastp_cmd} \
		--detect_adapter_for_pe --in1 ${fqs[0]} --in2 ${fqs[1]} \
		--out1 ${meta.id}_1${extension} --out2 ${meta.id}_2${extension} \
	rm ${fqs}
	"""
	else if (meta.single_end)
	"""
	${fastp_cmd} \
		--in1 ${fqs} --out1 ${meta.id}${extension} \
	rm ${fqs}
	"""
}

process priceseqfilter {
	label 'price'
	label 'cpu_medium'
	cpus 8
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
	def price_cmd = """PriceSeqFilter -a 8 -rnf 90 -rqf 85 0.98 -log c """
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
	fastqs = channel.fromPath(params.fastq_path)
		.map { fq -> [[id: fq.baseName() - ~/_[12]/], fq] }
	
	filter_barcodes(fastqs) | fastp | priceseqfilter
}
