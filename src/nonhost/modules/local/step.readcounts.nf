#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.fastqPath = "fastq/*/"
params.publishDir = "$PWD"
params.publishIntermediate = true
params.refIndexes = null //"Danio_rerio.GRCz11.108"
params.cleanupScript = ""
params.tmp = null
params.backupTmp = null
params.backupScratchHack = false
params.nxfUnstageHack = false

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
	tuple val(meta), path("counts/Aligned.out.bam"), emit: bam

	script:
	def STAR_READCOUNTS_CMD = """STAR --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
		--outFilterMismatchNoverLmax 0.04 \
		--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
		--alignIntronMin 20 --alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--outSAMtype BAM Unsorted --outSAMattributes NH HI NM MD \
		--genomeLoad NoSharedMemory --outReadsUnmapped None \
		--runThreadN ${task.cpus} """
	if (!meta.single_end)
	"""
	trap 'echo "Interrupt by external (OOM?), exiting."; exit 130' SIGINT

	${STAR_READCOUNTS_CMD} \
		--genomeDir ${indexes_dir2} \
		--readFilesCommand ${task.ext.gzipCmd} -dc \
		--outFileNamePrefix counts/ \
		--readFilesIn ${fqgz[0]} ${fqgz[1]}
	"""
	else if (meta.single_end)
	"""
	trap 'echo "Interrupt by external (OOM?), exiting."; exit 130' SIGINT

	${STAR_READCOUNTS_CMD} \
		--genomeDir ${indexes_dir2} \
		--readFilesCommand ${task.ext.gzipCmd} -dc \
		--outFileNamePrefix counts/ \
		--readFilesIn ${fqgz}
	"""
}


process sort_bam {
	label 'samtools'

	input:
	tuple val(meta), path(bam)

	output:
	tuple val(meta), path("Aligned.out.namesorted.bam"), emit: bam

	script:
	def mem = "${(task.memory.toMega() * 0.75).toLong().intdiv(task.cpus)}M"
	"""
	trap 'echo "Interrupt by external (OOM?), exiting."; exit 130' SIGINT

	samtools sort -m $mem -n -@ $task.cpus \
		-o Aligned.out.namesorted.bam $bam
	"""
}

process htseq_count {
	label 'htseq'

	input:
	tuple val(meta), path(sorted_bam)
	path gtf_noERCC // "Danio_rerio.GRCz11.108.gtf"

	output:
	tuple val(meta), path("counts.txt"), emit: counts
	tuple val(meta), path("htseq-mapping-extra.txt"), emit: summary
	tuple val(meta), path("counts-row.txt"), emit: countsRow

	script:
	"""
	htseq-count -r name -s no -f bam -m intersection-nonempty \
		$sorted_bam $gtf_noERCC > counts.txt
	rm $sorted_bam
	grep "^__" counts.txt > htseq-mapping-extra.txt
	grep -v "^__" counts.txt > counts.txt

	printf "Run," >counts-row-staging.txt
	cat counts.txt | cut -f1 | tr '\n' ',' | sed 's/,$//' >> counts-row-staging.txt
	printf "\n" >> counts-row-staging.txt
	printf "${meta.id}," >>counts-row-staging.txt
	cat counts.txt | cut -f2 | tr '\n' ',' | sed 's/,$//' >> counts-row-staging.txt
	printf "\n" >> counts-row-staging.txt
	mv counts-row-staging.txt counts-row.txt

	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
}

process feature_count {
	label 'featurecounts'

	input:
	tuple val(meta), path("${meta.id}.bam")
	path gtf_noERCC // "Danio_rerio.GRCz11.108.gtf"

	output:
	tuple val(meta), path("counts.txt"), emit: counts
	tuple val(meta), path("feature-counts.txt.summary"), emit: summary
	tuple val(meta), path("counts-row.txt"), emit: countsRow

	script:
	def p = "${meta.single_end ? '' : '-p --countReadPairs'}"
	"""
	sorted_bam="${meta.id}.bam"
	featureCounts $p -T ${task.cpus} -a $gtf_noERCC -o feature-counts.staging.txt \
		\$sorted_bam
	<feature-counts.staging.txt tail -n +2 | cut -f1,7 >counts.txt

	printf "Run," >counts-row-staging.txt
	<counts.txt tail -n +2 | cut -f1 | tr '\n' ',' | sed 's/,$//' >> counts-row-staging.txt
	printf "\n" >> counts-row-staging.txt
	printf "${meta.id}," >>counts-row-staging.txt
	<counts.txt tail -n +2 | cut -f2 | tr '\n' ',' | sed 's/,$//' >> counts-row-staging.txt
	printf "\n" >> counts-row-staging.txt
	mv counts-row-staging.txt counts-row.txt
	
	mv feature-counts.staging.txt.summary feature-counts.txt.summary
	rm feature-counts.staging.txt
	
	cleanup="${meta.cleanup}"
	${params.cleanupScript}
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
