#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def container_usage() {
	return """Containerization options:
-with-docker              use docker containers to run commands when possible
-with-singularity         use singularity containers to run commands when possible
                            (avoid most pre-installation procedure)
"""
}
def publish_usage() {
	return """Intermediate output options:
--publish_intermediate    if true, publish some intermediate step output
                            (default: false)
--publish_dir path        path to write useful intermediate output of each step
                            (default: \$PWD)
--publish_fastqs          whether to publish output from fastq dump
                            (default: true)
--publish_pricefiltered   whether to publish output from PRICE filtering
                            (default: true)
--publish_star            whether to publish output from STAR mapping
                            (default: true)
--publish_readcounts      whether to publish output from the read count pipeline
                            (default: true)
"""
}

def helpMessage() {
	log.info """
--accessions_list path    file containing sra accessions to process, one per line
                            (required)
--parallel_prefetch n     maximum sra prefetch downloads to run at once
                            (default: 10)
--help, -h                print this text and exit

${star_usage()}
${container_usage()}
${publish_usage()}
"""
}

params.help = false
params.h = false
params.accessions_list = ""

params.genome_size = null // must be filled
params.star_ref_indexes_ercc = null // "Danio_rerio.GRCz11.108.ERCC"
params.star_ref_indexes = null // "Danio_rerio.GRCz11.108"
params.hisat2_ref_indexes = null
// https://ftp.ensembl.org/pub/release-108/fasta/danio_rerio/dna/
params.ref_genome = "Danio_rerio.GRCz11.dna_sm.primary_assembly.fa"
// https://ftp.ensembl.org/pub/release-108/gtf/danio_rerio/
params.ref_genome_gtf = "Danio_rerio.GRCz11.108.gtf"
// https://tools.thermofisher.com/content/sfs/manuals/ERCC92.zip
params.ercc_fa = "ERCC92.fa"
params.ercc_gtf = "ERCC92.gtf"

params.parallel_downloads = 100
params.skip_host_counts = false
params.skip_hisat2 = false
params.hisat2_use_transcript = true

params.publish_dir = "$PWD"
params.publish_intermediate = false
params.publish_fastqs = true
params.publish_pricefiltered = true
params.publish_readcounts = true
params.publish_hisat2 = true
params.publish_star = true

params.star_index_gen_options = ""
params.hisat2_index_gen_options = ""
params.sra_prefetch_options = ""
params.fastq_dump_options = ""
params.fastp_options = ""
params.price_options = ""
params.star_count_options = "" // host read counts pipeline
params.samtools_sort_options = ""
params.htseq_count_options = ""
params.hisat2_options = "" // nonhost pipeline
params.star_options = ""

include {
	star_generate_indexes;
	hisat2_generate_indexes;
} from './modules/step.0.generate_indexes.nf' params(
	publish_dir: params.publish_dir,
	star_index_gen_options: params.star_index_gen_options,
	hisat2_index_gen_options: params.hisat2_index_gen_options
	hisat2_use_transcript: params.hisat2_use_transcript
)
include {
	prefetch;
	fastq_dump;
	filter_barcodes;
} from './modules/step.1.nf' params(
	parallel_downloads: params.parallel_downloads,
	publish_dir: params.publish_dir,
	publish_intermediate: params.publish_intermediate && params.publish_fastqs,
	sra_prefetch_options: params.sra_prefetch_options,
	fastq_dump_options: params.fastq_dump_options
)
include {
	fastp;
	priceseqfilter;
} from './modules/step.2.nf' params(
	publish_dir: params.publish_dir,
	publish_intermediate: params.publish_intermediate && params.publish_pricefiltered,
	fastp_options: params.fastp_options,
	price_options: params.price_options
)
include {
	star_counts;
	sort_bam;
	htseq_count;
} from './modules/step.readcounts.nf' params(
	genome_size: params.genome_size,
	publish_dir: params.publish_dir,
	publish_intermediate: params.publish_intermediate && params.publish_readcounts,
	star_count_options: params.star_count_options,
	samtools_sort_options: params.samtools_sort_options,
	htseq_count_options: params.htseq_count_options
)
include {
	hisat2;
	ensure_hisat2_indexes;
} from './modules/step.3.hisat2.nf' params(
	genome_size: params.genome_size,
	publish_dir: params.publish_dir,
	publish_intermediate: params.publish_intermediate && params.publish_hisat2,
	hisat2_options: params.hisat2_options
)
include {
	star_usage;
	star;
	ensure_star_indexes;
} from './modules/step.4.star.nf' params(
	genome_size: params.genome_size,
	publish_dir: params.publish_dir,
	publish_intermediate: params.publish_intermediate && params.publish_star,
	star_options: params.star_options
)


workflow {
	if (params.help || params.h || !params.accessions_list) {
		helpMessage()
		exit params.help || params.h ? 0 : 1
	}
	
	main:
	// step 0: generating hisat2 indexes
	hisat2_indexes = ensure_hisat2_indexes(params.hisat2_ref_indexes,
										   params.ref_genome,
										   params.ref_genome_gtf,
										   params.ercc_fa,
										   params.ercc_gtf)
	// step 0: generating STAR indexes
	(star_indexes, star_indexes2) = ensure_star_indexes(params.star_ref_indexes,
														params.star_ref_indexes_ercc,
														params.ref_genome,
														params.ref_genome_gtf,
														params.ercc_fa,
														params.ercc_gtf)


	// step 1: download and convert to fastq
	accessions = channel.fromPath(params.accessions_list).splitText()
		.map { acc -> [[id: acc.trim()], acc.trim()] }

	sra = prefetch(accessions)
		.map { meta, sra, reads, sra_size ->
			def new_meta = [id: meta.id,
							reads: reads.toInteger(),
							sra_size: sra_size.toInteger() ]
			[ new_meta, sra ]
		}
		
	fastqs = fastq_dump(sra)

	fastqs_2 = filter_barcodes(fastqs)
		.map { meta, fastq, median, count, fsize ->
			def new_meta = meta.clone()
			new_meta.reads = count.toInteger()
			new_meta.readlen = median.toInteger()
			new_meta.fastq_size = fsize.toInteger()
			new_meta.single_end = fastq.size() != 2
			[ new_meta, fastq ]
		}

	// step 2: adapter trimming & filtering
	filtered_fastqs = fastqs_2 | fastp | priceseqfilter

	// read counts path
	if (!params.skip_host_counts) {
		bam = star_counts(filtered_fastqs, star_indexes2)
		bam_sorted = sort_bam(bam)
		count = htseq_count(bam_sorted, file(params.ref_genome_gtf))
	}

	// step 3: hisat2
	if (!params.skip_hisat2) {
		unmapped_reads_1 = hisat2(filtered_fastqs, hisat2_indexes)
	} else {
		unmapped_reads_1 = filtered_fastqs
	}
	
	// step 4: STAR
	unmapped_reads_2 = star(unmapped_reads_1, star_indexes)

}
