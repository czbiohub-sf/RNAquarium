#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.fastq_path = "fastq/*/"
params.publish_dir = "$PWD"
params.publish_intermediate = true

process filter_barcodes {
	label 'median'

	input:
	path fastq

	output:
	path 'fastq/*.fastq'

	shell:
	'''
	acc=$(ls -AU fastq | head -1 | grep -Eo "[SED]RR[0-9]*")
	# unzipped may be cached
	if [ $(ls fastq/*.gz | wc -w) -gt 0 ]; then gunzip -q fastq/*.gz; fi

	# one set of reads, or two?
	if [[ $(ls -1 fastq | wc -l) -lt 2 ]]
	then # single
		#echo $(fastq-lengths median 1 fastq/${acc}.fastq)
		echo $acc SE
	else # possibly paired, but may be scRNAseq barcodes
		if [ $(fastq-lengths median fastq/${acc}_1.fastq) -lt 32 ] && \
			[ $(fastq-lengths median fastq/${acc}_2.fastq) -gt 80 ]
		then
			echo $acc scRNAseq
			rm fastq/${acc}_1.fastq # discard read 1 (cell barcode)
			mv fastq/${acc}_2.fastq fastq/${acc}.fastq
		else
			echo $acc PE
		fi
	fi
	'''
}

process fastp {
	label 'fastp'

	input:
	path fqs

	output:
	path '*.trimmed.fastq'

	script:
	def extension = ".trimmed.fastq"
	if (fqs.getClass() == nextflow.util.BlankSeparatedList)
	"""
	fastp --disable_quality_filtering --disable_length_filtering \
		--compression 6 --thread ${task.cpus} \
		--detect_adapter_for_pe --in1 ${fqs[0]} --in2 ${fqs[1]} \
		--out1 ${fqs[0].baseName}${extension} --out2 ${fqs[1].baseName}${extension} \
		--json fastp.json --html fastp.html
	rm ${fqs}
	"""
	else if (fqs.getClass() == nextflow.processor.TaskPath)
	"""
	fastp --disable_quality_filtering --disable_length_filtering \
		--compression 6 --thread ${task.cpus} \
		--in1 ${fqs} --out1 ${fqs.baseName}${extension} \
		--json fastp.json --html fastp.html
	rm ${fqs}
	"""
}

process priceseqfilter {
	label 'price'
	publishDir params.publish_dir, enabled: params.publish_intermediate
	
	input:
	path fqs

	output:
	path '*.PRICEfiltered.fastq.gz'
	
// -a thread, -log c concise output, -fp input, -op output
// -rnf 90:    90 percentage of nucleotides in a read that must be called
// -rqf 85 0.98:       85% of the read length must be 98% accurate (parameterized after CZID, Joe uses 95% of the read length 98% accurate)
	script:
	def extension = ".PRICEfiltered.fastq"
	if (fqs.getClass() == nextflow.util.BlankSeparatedList)
	"""
	PriceSeqFilter -a ${task.cpus} -rnf 90 -rqf 85 0.98 -log c \
		-fp ${fqs} \
		-op ${fqs[0].baseName}${extension} ${fqs[1].baseName}${extension}
	gzip *${extension}
	rm ${fqs}
	"""
	else if (fqs.getClass() == nextflow.processor.TaskPath)
	"""
	PriceSeqFilter -a ${task.cpus} -rnf 90 -rqf 85 0.98 -log c \
		-f ${fqs} -o ${fqs.baseName}${extension}
	gzip *${extension}
	rm ${fqs}
	"""
}

workflow {
	fastqs = channel.fromPath(params.fastq_path)
	filter_barcodes(fastqs) | fastp | priceseqfilter
}
