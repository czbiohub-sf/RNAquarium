#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.fastqPath = "fastq/*/"
params.publishDir = "$PWD"
params.publishIntermediate = true
params.refIndexes = null //"Danio_rerio.GRCz11.108"

params.starCountOptions = ""
params.samtoolsSortOptions = ""
params.htseqCountOptions = ""

params.metaIn = 'step_2_sheet.csv'
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
		$params.starCountOptions"""
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
	samtools sort -m $mem -n -@ $task.cpus $params.samtoolsSortOptions \
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
		$params.htseqCountOptions $sorted_bam $gtf_noERCC > htseq-count.txt
	rm $sorted_bam
	"""
}

workflow {
	indexes2 = params.refIndexes
	
	fastqs = channel.fromFilePairs("$params.fastqPath/[SED]RR*_?[12]?.fastq",
								   size: -1)
	
	bam = star_counts(fastqs, indexes2)
	bam_sorted = sort_bam(bam)
	count = htseq_count(bam_sorted, params.gtfNoErcc)
}