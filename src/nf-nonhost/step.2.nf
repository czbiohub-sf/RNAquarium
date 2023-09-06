#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.fastq_path = "fastq/*/"

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

/*	stub:
	'''
	find fastq -name "*.gz" -exec sh -c 'mv "$1" "${1%.gz}"' _ {} \\;
	'''*/
}


process fastp_pe {
	debug true
	label 'fastp'
	// cpus 4
	// memory '24GB'
	// todo: scratch directive for hpc

	input:
	tuple path(fq1), path(fq2)

	output:
	path '*.trimmed.fastq'
	
	script:
	def extension = ".trimmed.fastq"
	"""
	fastp --disable_quality_filtering --disable_length_filtering \
		--compression 6 --thread ${task.cpus} --detect_adapter_for_pe \
		--in1 ${fq1} --in2 ${fq2} \
		--out1 ${fq1.baseName}${extension} --out2 ${fq2.baseName}${extension} \
		--json fastp.json --html fastp.html
	rm ${fq1} ${fq2}
	"""
}
process fastp_se {
	label 'fastp'

	input:
	path fq

	output:
	path '*.trimmed.fastq'

	script:
	def extension = ".trimmed.fastq"
	"""
	fastp --disable_quality_filtering --disable_length_filtering \
		--compression 6 --thread ${task.cpus} \
		--in1 ${fq} --out1 ${fq.baseName}${extension} \
		--json fastp.json --html fastp.html
	rm ${fq}
	"""
}

workflow {
	fastqs = channel.fromPath(params.fastq_path)
	fastqs.view { "$it" }
}
