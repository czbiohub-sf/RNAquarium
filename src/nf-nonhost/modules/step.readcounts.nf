#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.fastq_path = "fastq/*/"
params.publish_dir = "$PWD"
params.publish_intermediate = true
params.ref_indexes = null //"Danio_rerio.GRCz11.108"

params.star_count_options = ""
params.samtools_sort_options = ""
params.htseq_count_options = ""

params.meta_in = 'step_2_sheet.csv'
include {
	LOAD_METASHEET;
} from './utils.nf'

process star_counts {
	label 'star'
	
	input:
	tuple val(meta), path(fqgz)
	path indexes_dir2 // "Danio_rerio.GRCz11.108"
	
	output:
	tuple val(meta), path("counts/Aligned.out.bam")

	script:
	def STAR_READCOUNTS_CMD = """STAR --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
		--outFilterMismatchNoverLmax 0.04 \
		--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
		--alignIntronMin 20 --alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--outSAMtype BAM Unsorted --outSAMattributes NH HI NM MD \
		--genomeLoad NoSharedMemory --outReadsUnmapped None \
		--runThreadN ${task.cpus} \
		$params.star_count_options"""
	if (!meta.single_end)
	"""
	${STAR_READCOUNTS_CMD} \
		--genomeDir ${indexes_dir2} \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix counts/ \
		--readFilesIn ${fqgz[0]} ${fqgz[1]}
	"""
	else if (meta.single_end)
	"""
	${STAR_READCOUNTS_CMD} \
		--genomeDir ${indexes_dir2} \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix counts/ \
		--readFilesIn ${fqgz}
	"""
}


process sort_bam {
	label 'samtools'

	input:
	tuple val(meta), path(bam)

	output:
	tuple val(meta), path("Aligned.out.namesorted.bam")

	script:
	def mem = "${task.memory.toGiga()}G"
	"""
	samtools sort -m $mem -n -@ $task.cpus $params.samtools_sort_options \
		-o Aligned.out.namesorted.bam $bam
	"""
}

process htseq_count {
	label 'htseq'
	
	input:
	tuple val(meta), path(sorted_bam)
	path gtf_noERCC // "Danio_rerio.GRCz11.108.gtf"

	output:
	tuple val(meta), path("htseq-count.txt")

	script:
	"""
	htseq-count -r name -s no -f bam -m intersection-nonempty \
		$params.htseq_count_options $sorted_bam $gtf_noERCC > htseq-count.txt
	rm $sorted_bam
	"""
}

workflow {
	indexes2 = params.ref_indexes
	
	fastqs = channel.fromFilePairs("$params.fastq_path/[SED]RR*_?[12]?.fastq",
								   size: -1)
	
	bam = star_counts(fastqs, indexes2)
	bam_sorted = sort_bam(bam)
	count = htseq_count(bam_sorted, params.gtf_no_ercc)
}