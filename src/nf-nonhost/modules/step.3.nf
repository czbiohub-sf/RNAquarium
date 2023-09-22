#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.fastq_path = "fastq/*/"
params.publish_dir = "$PWD"
params.publish_intermediate = true
params.ref_indexes_ercc = null // "Danio_rerio.GRCz11.108.ERCC"
params.ref_indexes = null //"Danio_rerio.GRCz11.108"

// https://ftp.ensembl.org/pub/release-108/fasta/danio_rerio/dna/
// https://ftp.ensembl.org/pub/release-108/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz
params.ref_genome = "Danio_rerio.GRCz11.dna_sm.primary_assembly.fa"
// https://ftp.ensembl.org/pub/release-108/gtf/danio_rerio/
// https://ftp.ensembl.org/pub/release-108/gtf/danio_rerio/Danio_rerio.GRCz11.108.gtf.gz
params.gtf_no_ercc = "Danio_rerio.GRCz11.108.gtf"
// https://tools.thermofisher.com/content/sfs/manuals/ERCC92.zip
params.ercc = "ERCC92.fa"
params.ercc_gtf = "ERCC92.gtf"

process generate_indexes {
	label 'star'
	cpus 16
	cache true
	publishDir params.publish_dir

	input:
	path ref_genome_fa
	path ref_genome_gtf
	path ERCC_fa
	path ERCC_gtf

	output:
	path("indexes.ERCC"), emit: ercc_indexes
	path("indexes"), emit: indexes

	script:
	"""
	cat $ref_genome_fa $ERCC_fa > dna_sm.primary_assembly_ERCC.fa
	cat $ref_genome_gtf $ERCC_gtf > indexes_ERCC.gtf

	STAR --runMode genomeGenerate --runThreadN $task.cpus --genomeFastaFiles \
		$ref_genome_fa --sjdbGTFfile $ref_genome_gtf \
		--genomeDir	./indexes

	STAR --runMode genomeGenerate --runThreadN $task.cpus --genomeFastaFiles \
		dna_sm.primary_assembly_ERCC.fa --sjdbGTFfile indexes_ERCC.gtf \
		--genomeDir	./indexes.ERCC
	"""
}
	
process star {
	label 'star'
	cpus 16
	publishDir "$params.publish_dir/STAR_out/$meta.id", enabled: params.publish_intermediate

	input:
	tuple val(meta), path(fqgz)
	path indexes_dir // "Danio_rerio.GRCz11.108.ERCC"

	output:
	path "?E/Unmapped.out.mate?.gz"

	script:
	def star_cmd = """STAR --outFilterMultimapNmax 99999 --outFilterMismatchNmax 999 \
		--outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5 \
		--outSAMmode None --quantMode GeneCounts \
		--clip3pNbases 0 --limitOutSJcollapsed 200000000 \
		--genomeLoad NoSharedMemory --outReadsUnmapped Fastx \
		--runThreadN ${task.cpus}"""
	if (!meta.single_end)
	"""
	${star_cmd} \
		--genomeDir ${indexes_dir} \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix PE/ \
		--readFilesIn ${fqgz[0]} ${fqgz[1]}

	gzip Unmapped.out.mate*
	"""
	else if (meta.single_end)
	"""
	${star_cmd} \
		--genomeDir ${indexes_dir} \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix SE/ \
		--readFilesIn ${fqgz}

	gzip Unmapped.out.mate1
	"""
}

process star_counts {
	label 'star'
	cpus 16
	publishDir "$params.publish_dir/STAR_out/counts/$meta.id", enabled: params.publish_intermediate

	input:
	tuple val(meta), path(fqgz)
	path indexes_dir2 // "Danio_rerio.GRCz11.108"
	
	output:
	tuple val(meta), path("counts/Aligned.out.bam")

	script:
	def star_readcounts_cmd = """STAR --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
		--outFilterMismatchNoverLmax 0.04 \
		--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
		--alignIntronMin 20 --alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--outSAMtype BAM Unsorted --outSAMattributes NH HI NM MD \
		--genomeLoad NoSharedMemory --outReadsUnmapped None \
		--runThreadN ${task.cpus}"""
	if (!meta.single_end)
	"""
	${star_readcounts_cmd} \
		--genomeDir ${indexes_dir2} \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix counts/ \
		--readFilesIn ${fqgz[0]} ${fqgz[1]}
	"""
	else if (meta.single_end)
	"""
	${star_readcounts_cmd} \
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
	"""
	samtools sort -m 45G -n -o -@ $task.cpus -o Aligned.out.namesorted.bam $bam
	"""
}

process htseq_count {
	label 'htseq'
	publishDir "$params.publish_dir/counts/$meta.id", enabled: params.publish_intermediate

	input:
	tuple val(meta), path(sorted_bam)
	path gtf_noERCC // "Danio_rerio.GRCz11.108.gtf"

	output:
	tuple val(meta), path("htseq-count.txt")

	script:
	"""
	htseq-count -r name -s no -f bam -m intersection-nonempty \
		$sorted_bam	$gtf_noERCC > htseq-count.txt
	rm sorted_bam
	"""
}

workflow {
	// TODO: name better
	if (params.ref_indexes && params.ref_indexes_ercc) {
		indexes = params.ref_indexes_ercc
		indexes2 = params.ref_indexes
	} else {
		(indexes, indexes2) = generate_indexes(file(params.ref_genome),
											   file(params.gtf_no_ercc),
											   file(params.ercc),
											   file(params.indexes_ercc))
	}

	fastqs = channel.fromFilePairs("$params.fastq_path/[SED]RR*_?[12]?.fastq",
								   size: -1)
	
	unmapped_reads = star(fastqs, indexes)
	bam = star_counts(fastqs, indexes2)
	bam_sorted = sort_bam(bam)
	count = htseq_count(bam_sorted, params.gtf_no_ercc)
}
