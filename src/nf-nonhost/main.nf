#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def container_usage() {
	return """Containerization options:
-profile docker           use docker containers to run commands when possible
-profile singularity      use singularity containers to run commands when possible
							(avoid most pre-installation procedure)
"""
}

params.accessionList = ""
params.fastqPath = null // "$PWD/fastq"
params.parallelDownloads = 100
params.skipHostCounts = false
params.skipHisat = false
params.hisatUseTranscript = true
params.help = false

params.genomeSize = null // must be filled
params.starRefIndexesErcc = null // "Danio_rerio.GRCz11.108.ERCC"
params.starRefIndexes = null // "Danio_rerio.GRCz11.108"
params.hisatRefIndexes = null
params.bowtieRefIndexes = null
params.gsnapRefIndexes = null
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
params.publishBowtie = true
params.publishDedup = true
params.publishGsnap = true

params.starUseSharedMem = false
params.starThreadsSmall = 4
params.starThreadsLarge = 16
params.retainMixed = true
params.dedupPercentLen = 100
params.dedupMinLen = 20

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'
include { } from './modules/local/step.0.generate_indexes.nf' params(
	publishDir: params.publishDir,
	hisatUseTranscript: params.hisatUseTranscript
)
include {
	prefetch;
	fastq_dump;
	check_direct_fastqs;
	filter_barcodes;
} from './modules/local/step.1.nf' params(
	parallelDownloads: params.parallelDownloads,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishFastqs,
)
include {
	fastp;
	priceseqfilter;
} from './modules/local/step.2.nf' params(
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishPricefiltered,
)
include {
	star_counts;
	sort_bam;
	htseq_count;
} from './modules/local/step.readcounts.nf' params(
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishReadcounts,
)
include {
	hisat2;
	ensure_hisat2_indexes;
} from './modules/local/step.3.hisat2.nf' params(
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishHisat,
)
include {
	star;
	ensure_star_indexes;
} from './modules/local/step.4.star.nf' params(
	genomeSize: params.genomeSize,
	starUseSharedMem: params.starUseSharedMem,
	starThreadsSmall: params.starThreadsSmall,
	starThreadsLarge: params.starThreadsLarge,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishStar,
)
include {
	bowtie2;
	process_bowtie2_sam;
	bowtie2_filter_by_names;
	ensure_bowtie2_indexes;
} from './modules/local/step.5.bowtie.nf' params(
	retainMixed: params.retainMixed, // retain mixed (xor) mate align cases
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishBowtie,
)
include {
	dedup;
} from './modules/local/step.6.dedup.nf' params(
	percentLen: params.dedupPercentLen,
	minLen: params.dedupMinLen,
	pubishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishDedup,
)
include {
	gsnap;
	process_gsnap_sam;
	gsnap_filter_by_names;
	gsnap_skip;
	ensure_gsnap_indexes;
} from './modules/local/step.7.gsnap.nf' params(
	retainMixed: params.retainMixed, // retain mixed (xor) mate align cases
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishGsnap,
)
include {
	stats_csv;
} from './modules/local/stats.nf' params(
	skipHisat: params.skipHisat,
)


workflow {
	// parameter validation:
	// (in addition to format validation schema.json)
	// at least one of accessionList or fastqPath MUST be specified
	// genomeSize MUST be specified
	// parallelDownloads MUST be >= 1
	// refGenomeGtf MUST be specified if skipHostCounts is NOT specified
	// starRefIdx, starRefIdxErcc, MUST be present IF:
	//    refGenome, refGenomeGtf, erccFa, OR erccGtf are not
	// if hisatUse  hisatRefIdx
	if (params.help || !(params.accessionList || params.fastqPath)) {
		log.info paramsHelp("""PATH=\$PATH:\$PWD/bin nextflow run main.nf
	--accession-list sras.txt --ref-genome Danio_rerio.GRCz11.dna_sm.primary_assembly.fa
	--ref-genome-gtf Danio_rerio.GRCz11.108.gtf --ercc-fa ERCC92.fa --ercc-gtf ERCC92.gtf --genome-size 1396431182 -profile slurm,apptainer
	--tmp=/tmp/""")
		exit params.help ? 0 : 1
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
	validateParameters()
	StringBuilder param_info = new StringBuilder()
	for (e in params) {
		if (e.key.equals(e.key.toLowerCase()) && !(e.key in ['help', 'h']))
			param_info.append("${e.key.padLeft(23)}:\t$e.value\n")
	}
	//log.info param_info.toString()
	log.info paramsSummaryLog(workflow)
	if (!params.tmp) log.warn "--tmp=<path> not specified. using /tmp/\n(choose a scratch space appropriate for many very large files)"

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
	// step 0: generating bowtie2 indexes
	bowtie2_indexes = ensure_bowtie2_indexes(params.bowtieRefIndexes,
											 params.refGenome,
											 params.erccFa)
	// step 0: generating gsnap indexes
	gsnap_indexes = ensure_gsnap_indexes(params.gsnapRefIndexes,
										 params.refGenome,
										 params.erccFa)

	// step 1: download and convert to fastq
	// find existing fastqs
	accessions = Channel.fromPath(params.accessionList, type: 'file').splitText()
		.map { acc -> acc.trim() }

	if (params.fastqPath && file(params.fastqPath).exists()) {
		try {
			direct_fastqs = Channel.fromPath("$params.fastqPath/*", type: 'dir')
				.map { path -> // need to think about this more, failure handling?
					def new_meta = [ id: path.getSimpleName(),
									sra_size: files("$path/*.fastq")[0].size() ]
					[ new_meta, path ]
				}
			direct_fastq_ids = direct_fastqs
				.map { meta, _ ->
					[ meta.id, true ]
				}
			// remove ids that exist in pre-dumped fastq path from accessions list
			accessions = accessions
				.join(direct_fastq_ids, remainder: true, by: 0)
				.filter { key, v2 -> !v2 }
				.map { key, _ ->
					[ [id: key], key ]
				}
		} catch (Exception e) {
			log.error "--fastq-path $params.fastqPath is not folders of fastq?\n$e"
			exit(1)
		}
	} else {
		direct_fastqs = Channel.empty()
		accessions = accessions
			.map { key ->
				[ [id: key], key ]
			}
	}

	// prefetch SRAs by remaining accessions
	sra = prefetch(accessions).sra
		.map { meta, sra, reads, sra_size ->
			def new_meta = [id: meta.id,
							reads: reads.toLong(),
							sra_size: sra_size.toLong() ]
			[ new_meta, sra ]
	}

	// and convert SRA to fastq
	fastqs = fastq_dump(sra).mates
		.mix(check_direct_fastqs(direct_fastqs)) // , merging any existing fastqs

	// heuristic filter scRNAseq barcode files and add metadata for resource optimization
	filter_barcodes(fastqs)
		.map { meta, fastq, median, count, fsize ->
			def new_meta = meta.clone()
			new_meta.reads = count.toLong()
			new_meta.readlen = median.toLong()
			new_meta.fastq_size = fsize.toLong()
			new_meta.single_end = fastq.size() != 2
			[ new_meta, fastq ]
		}
		.branch { // empty/insignificant runs (by trimming, qc, or host mapping) should drop out)
			ok: { meta, fastq -> 
				meta.single_end ? fastq.size() > 132 : fastq[0].size() > 132
				}
			dropouts: true
		}
	    .set { filter_barcodes_result }

	// step 2: adapter trimming & filtering
	fastp(filter_barcodes_result.ok)
	fastp.out.mates
		.branch { // empty/insignificant runs (by trimming, qc, or host mapping) should drop out)
			ok: { meta, fastq -> 
				meta.single_end ? fastq.size() > 132 : fastq[0].size() > 132
				}
			dropouts: true
		}
		.set { fastp_result }

	priceseqfilter(fastp_result.ok).mates
		.branch { // empty/insignificant runs (by trimming, qc, or host mapping) should drop out)
			ok: { meta, fastq -> 
				meta.single_end ? fastq.size() > 132 : fastq[0].size() > 132
			}
			dropouts: true
		}
		.set { priceseqfilter_result }

	// host read counts path
	if (!params.skipHostCounts) {
		star_counts(priceseqfilter_result.ok, star_indexes2)
		sort_bam(star_counts.out.bam)
		count = htseq_count(sort_bam.out.bam, file(params.refGenomeGtf))
	}

	// step 3: hisat2
	if (!params.skipHisat) {
		hisat2(priceseqfilter_result.ok, hisat2_indexes)
		unmapped_reads_1 = hisat2.out.mates
	} else {
		unmapped_reads_1 = filtered_fastqs
	}

	// step 4: STAR
	star(unmapped_reads_1, star_indexes)

	// step 5: bowtie2
	bowtie2(star.out.mates, bowtie2_indexes)
	process_bowtie2_sam(bowtie2.out.sam)
	bowtie2_filter_by_names(process_bowtie2_sam.out.names.join(star.out.mates))

	// step 6: deduplication
	bowtie2_filter_by_names.out.mates
		.branch { // empty/insignificant runs (by trimming, qc, or host mapping) should drop out)
			ok: { meta, fastq -> 
				meta.single_end ? fastq.size() > 132 : fastq[0].size() > 132
			}
			dropouts: true
		}
		.set { bowtie2_result }

	dedup(bowtie2_result.ok.join(star.out.stats))

	// step 7: gsnap
	gsnap(dedup.out.mates, gsnap_indexes)
	gsnap.out.sam
		.branch {
			ok: { meta, sam -> sam.size() > 0 }
			fails: true
		}
	    .set { gsnap_result }
	process_gsnap_sam(gsnap_result.ok)
	gsnap_filter_by_names(process_gsnap_sam.out.names.join(dedup.out.mates))

	// if gsnap fails a run for any non-oom reason, keep the nonhost reads from before that.
	gsnap_skip(gsnap_result.fails.join(dedup.out.mates))


	// rekey in preparation for join
	// questionable to use prefetch / fastq_dump stats b/c not present for direct fastq
	// but we DO
	//sra_stats = sra.out.stats.map { meta, sra -> [ meta.id, sra ] }
	//fastq_stats = fastqs.out.stats.map { meta, fastq -> [ meta.id, sra ] }
	meta_stats    = filter_barcodes_result.ok.map { meta, fastq -> [ meta.id, meta ] }
	fastp_stats   = fastp.out.stats.map           { meta, stats -> [ meta.id, stats ] }
	price_stats   = priceseqfilter.out.stats.map  { meta, stats -> [ meta.id, stats ] }
	if (!params.skipHisat) {
		hisat2_stats  = hisat2.out.stats.map  { meta, stats -> [ meta.id, stats ] }
	} else {
		hisat2_stats = [ meta_stats[0], "n/a" ]
	}
	star_stats    = star.out.stats.map    { meta, stats -> [ meta.id, stats ] }
	bowtie2_stats = process_bowtie2_sam.out.stats.map { meta, stats -> [ meta.id, stats ] }
	dedup_stats   = dedup.out.stats.map   { meta, stats -> [ meta.id, stats ] }
	gsnap_stats   = process_gsnap_sam.out.stats.map { meta, stats -> [ meta.id, stats ] }

	all_stats = meta_stats.join(fastp_stats)
			  .join(price_stats)
			  .join(hisat2_stats)
			  .join(star_stats)
			  .join(bowtie2_stats)
			  .join(dedup_stats)
			  .join(gsnap_stats)

	stats_csv(all_stats)
		.collectFile(name: "stats-${params.timestamp}.csv", keepHeader: true, skip: 1, storeDir: "${params.publishDir}/stats/")
}


