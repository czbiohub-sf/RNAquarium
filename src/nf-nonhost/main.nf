#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def container_usage() {
	return """Containerization options:
-profile docker           use docker containers to run commands when possible
-profile singularity      use singularity containers to run commands when possible
                            (avoid most pre-installation procedure)
"""
}
def publish_usage() {
	return """Intermediate output options:
--publish-intermediate    if true, publish some intermediate step output
                            (default: false)
--publish-dir path        path to write useful intermediate output of each step
                            (default: \$PWD)
--publish-fastqs          whether to publish output from fastq dump
                            (default: true)
--publish-pricefiltered   whether to publish output from PRICE filtering
                            (default: true)
--publish-star            whether to publish output from STAR mapping
                            (default: true)
--publish-readcounts      whether to publish output from the read count pipeline
                            (default: true)
"""
}

def helpMessage() {
	log.info """
--accession-list path    file containing sra accessions to process, one per line
                            (required)
--parallel-prefetch n     maximum sra prefetch downloads to run at once
                            (default: 100)
--help, -h                print this text and exit
--genome-size n           genome size (approximate), in bytes
${star_usage()}
${container_usage()}
${publish_usage()}
"""
}

params.accessionList = ""
params.fastqPath = "$PWD/fastq"
params.parallelDownloads = 100
params.skipHostCounts = false
params.skipHisat = false
params.hisatUseTranscript = false
params.help = false
params.h = false

params.genomeSize = null // must be filled
params.starRefIndexesErcc = null // "Danio_rerio.GRCz11.108.ERCC"
params.starRefIndexes = null // "Danio_rerio.GRCz11.108"
params.hisatRefIndexes = null
// https://ftp.ensembl.org/pub/release-108/fasta/danio_rerio/dna/
params.refGenome = "Danio_rerio.GRCz11.dna_sm.primary_assembly.fa"
// https://ftp.ensembl.org/pub/release-108/gtf/danio_rerio/
params.refGenomeGtf = "Danio_rerio.GRCz11.108.gtf"
// https://tools.thermofisher.com/content/sfs/manuals/ERCC92.zip
params.erccFa = "ERCC92.fa"
params.erccGtf = "ERCC92.gtf"

params.publishDir = "$PWD"
params.publishIntermediate = false
params.publishFastqs = true
params.publishPricefiltered = true
params.publishReadcounts = true
params.publishHisat = true
params.publishStar = true

params.starIndexGenOptions = ""
params.hisatIndexGenOptions = ""
params.sraPrefetchOptions = ""
params.fastqDumpOptions = ""
params.fastpOptions = ""
params.priceOptions = ""
params.starCountOptions = "" // host read counts pipeline
params.samtoolsSortOptions = ""
params.htseqCountOptions = ""
params.hisatOptions = "" // nonhost pipeline
params.starOptions = ""

include {
	star_generate_indexes;
	hisat2_generate_indexes;
} from './modules/step.0.generate_indexes.nf' params(
	publishDir: params.publishDir,
	starIndexGenOptions: params.starIndexGenOptions,
	hisatIndexGenOptions: params.hisatIndexGenOptions,
	hisatUseTranscript: params.hisatUseTranscript
)
include {
	prefetch;
	fastq_dump;
	filter_barcodes;
} from './modules/step.1.nf' params(
	parallelDownloads: params.parallelDownloads,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishFastqs,
	sraPrefetchOptions: params.sraPrefetchOptions,
	fastqDumpOptions: params.fastqDumpOptions
)
include {
	fastp;
	priceseqfilter;
} from './modules/step.2.nf' params(
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishPricefiltered,
	fastpOptions: params.fastpOptions,
	priceOptions: params.priceOptions
)
include {
	star_counts;
	sort_bam;
	htseq_count;
} from './modules/step.readcounts.nf' params(
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishReadcounts,
	starCountOptions: params.starCountOptions,
	samtoolsSortOptions: params.samtoolsSortOptions,
	htseqCountOptions: params.htseqCountOptions
)
include {
	hisat2;
	ensure_hisat2_indexes;
} from './modules/step.3.hisat2.nf' params(
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishHisat,
	hisatOptions: params.hisatOptions
)
include {
	star_usage;
	star;
	ensure_star_indexes;
} from './modules/step.4.star.nf' params(
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishStar,
	starOptions: params.starOptions
)


workflow {
	// parameter validation:
	// at least one of accessionList or fastqPath MUST be specified
	// genomeSize MUST be specified (todo?: if starRefIndexesErcc not present)
	// parallelDownloads MUST be >= 1
	// refGenomeGtf MUST be specified if skipHostCounts is NOT specified
	// starRefIdx, starRefIdxErcc MUST be present
	//    if refGenome, refGenomeGtf, erccFa, OR erccGtf are not
	// if hisatUse  hisatRefIdx 
	if (params.help || params.h || !(params.accessionList || params.fastqPath)) {
		helpMessage()
		exit params.help || params.h ? 0 : 1
	} else if (!params.genomeSize || params.genomeSize <= 0) {
		log.error "--genome-size must be specified (approximate, in bytes)"
		exit 1
	} else if (params.parallelDownloads <= 0) {
		log.error "--parallel-downloads must be >= 1"
		exit 1
	} else if (!params.skipHostCounts && !file(params.refGenomeGtf).exists()) {
		log.error "--ref-genome-gtf annotations are required for host read counts"
		exit 1
	}

	StringBuilder param_info = new StringBuilder()
	for (e in params) {
		if (e.key.equals(e.key.toLowerCase()) && !(e.key in ['help', 'h']))
			param_info.append("${e.key.padLeft(23)}:\t$e.value\n")
	}
	log.info param_info.toString()
	
	main:
	// step 0: generating hisat2 indexes
	hisat2_indexes = ensure_hisat2_indexes(params.hisatRefIndexes,
										   params.refGenome,
										   params.refGenomeGtf,
										   params.erccFa,
										   params.erccGtf)
	// step 0: generating STAR indexes
	(star_indexes, star_indexes2) = ensure_star_indexes(params.starRefIndexes,
														params.starRefIndexesErcc,
														params.refGenome,
														params.refGenomeGtf,
														params.erccFa,
														params.erccGtf)

	// step 1: download and convert to fastq
	// find existing fastqs
	if (params.fastqPath) {
		direct_fastqs = Channel.fromPath("$params.fastqPath/*", type: 'dir')
			.map { path ->
				def new_meta = [ id: path.getSimpleName(),
								 sra_size: files("$path/*.fastq")[0].size() ]
				[ new_meta, path ]
			}
		direct_fastq_ids = direct_fastqs
			.map { meta, _ ->
				[ meta.id, true ]
			}
	} else {
		direct_fastq_ids = Channel.empty()
		direct_fastqs = Channel.empty()
	}

	// remove ids that exist in pre-dumped fastq path from accessions list
	accessions = Channel.fromPath(params.accessionList, type: 'file').splitText()
		.map { acc -> acc.trim() }
		.join(direct_fastq_ids, remainder: true, by: 0)
		.filter { key, v2 -> !v2 }
		.map { key, _ ->
			[ [id: key], key ]
		}
	
	// prefetch SRAs by remaining accessions
	sra = prefetch(accessions)
		.map { meta, sra, reads, sra_size ->
			def new_meta = [id: meta.id,
							reads: reads.toLong(),
							sra_size: sra_size.toLong() ]
			[ new_meta, sra ]
 	}

	// and convert SRA to fastq
	fastqs = fastq_dump(sra)
		.mix(direct_fastqs) // , merging any existing fastqs

	// heuristic filter scRNAseq barcode files and add metadata for resource optimization
	fastqs_2 = filter_barcodes(fastqs)
		.map { meta, fastq, median, count, fsize ->
			def new_meta = meta.clone()
			new_meta.reads = count.toLong()
			new_meta.readlen = median.toLong()
			new_meta.fastq_size = fsize.toLong()
			new_meta.single_end = fastq.size() != 2
			[ new_meta, fastq ]
		}

	// step 2: adapter trimming & filtering
	filtered_fastqs = fastqs_2 | fastp | priceseqfilter

	// host read counts path
	if (!params.skipHostCounts) {
		bam = star_counts(filtered_fastqs, star_indexes2)
		bam_sorted = sort_bam(bam)
		count = htseq_count(bam_sorted, file(params.refGenomeGtf))
	}

	// step 3: hisat2
	if (!params.skipHisat) {
		unmapped_reads_1 = hisat2(filtered_fastqs, hisat2_indexes)
	} else {
		unmapped_reads_1 = filtered_fastqs
	}
	
	// step 4: STAR
	unmapped_reads_2 = star(unmapped_reads_1, star_indexes)

}
