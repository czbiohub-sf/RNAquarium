#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.help = false
params.h = false
params.accessions_list = ""

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

params.parallel_downloads = 10
params.publish_dir = "$PWD"
params.publish_intermediate = false
params.publish_fastqs = true
params.publish_pricefiltered = true
params.publish_star = true

include {
	prefetch;
	fastq_dump;
} from './modules/step.1.nf' params(
	parallel_downloads: params.parallel_downloads,
	publish_dir: params.publish_dir,
	publish_intermediate: params.publish_intermediate && params.publish_fastqs
)
include {
	filter_barcodes;
	fastp;
	priceseqfilter;
} from './modules/step.2.nf' params(
	publish_dir: params.publish_dir,
	publish_intermediate: params.publish_intermediate && params.publish_pricefiltered
)
include {
	generate_indexes;
	star;
	star_counts;
	sort_bam;
	htseq_count;
} from './modules/step.3.nf' params(
	publish_dir: params.publish_dir,
	publish_intermediate: params.publish_intermediate && params.publish_star
)



def helpMessage() {
	log.info """
	--accessions_list path    file containing sra accessions to process, one per line
	                            (required)
	--parallel_prefetch n     maximum sra prefetch downloads to run at once
	                            (default: 10)
	--publish_dir path        path to write useful intermediate output of each step to
	--help, -h                print this text and exit
	
	-with-docker              use docker containers to run commands
	-with-singularity         use singularity containers to run commands
                                (avoid most pre-installation procedure)
	"""
}

workflow {
	if (params.help || params.h || !params.accessions_list) {
		helpMessage()
		exit params.help || params.h ? 0 : 1
	}
	
	main:
	// step 1
	accessions = channel.fromPath(params.accessions_list).splitText()
		.map { acc -> [[id: acc.trim()], acc.trim()] }
	fastqs = accessions | prefetch | fastq_dump

	// step 2
	fastqs_2 = filter_barcodes(fastqs)
		.map { meta, fastq ->
			def fmeta = meta
			if (fastq.size() == 2) {
				fmeta.single_end = false
			} else {
				fmeta.single_end = true
			}
			[ fmeta, fastq ]
		}
	filtered_fastqs = fastqs_2 | fastp | priceseqfilter


	if (params.ref_indexes && params.ref_indexes_ercc) {
		indexes = params.ref_indexes_ercc
		indexes2 = params.ref_indexes
	} else {
		(indexes, indexes2) = generate_indexes(file(params.ref_genome),
											   file(params.gtf_no_ercc),
											   file(params.ercc),
											   file(params.ercc_gtf))
	}

	unmapped_reads = star(filtered_fastqs, indexes)
	bam = star_counts(filtered_fastqs, indexes2)
	bam_sorted = sort_bam(bam)
	count = htseq_count(bam_sorted, params.gtf_no_ercc)

}
