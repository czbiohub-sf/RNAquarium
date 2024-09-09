#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.accessionList = ""
params.fastqPath = null // "$PWD/fastq"
params.parallelDownloads = 100
params.skipHostCounts = false
params.skipHisat = false
params.help = false

params.genomeSize = null // must be filled
params.starRefIndexesErcc = null
params.starRefIndexes = null
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
params.publishQCfiltered = true
params.publishReadcounts = true
params.publishHisat = true
params.publishStar = true
params.publishBowtie = true
params.publishDedup = true
params.publishGsnap = true

params.tmp = null
params.backupTmp = null
params.backupScratchHack = false
params.nxfUnstageHack = false

params.extraAdapters = "$PWD/extra-adapters.fasta"
params.hisatUseTranscript = true
params.starSjdbOverhang = 100
params.starUseSharedMem = false
params.starThreadsSmall = 4
params.starThreadsLarge = 16
params.retainMixed = true
params.dedupPercentLen = 100
params.dedupMinLen = 20

// https://pirl.unc.edu/blog/tricking-nextflows-caching-system-to-drastically-reduce-storage-usage
params.cleanupIntermediate = true
cleanupScript = params.cleanupIntermediate ? """for file in \$cleanup; do
if [ -e \$file ]; then
	size=`stat --printf="%s" \$file`
	atime=`stat --printf="%X" \$file`
	mtime=`stat --printf="%Y" \$file`

	# Make the file size 0 and set as a sparse file
	> \$file
	truncate -s \$size \$file
	# Reset the timestamps on the file
	touch -a -d @\$atime \$file
	touch -m -d @\$mtime \$file
fi
done""" : ""

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'
include { cleanup_branched; } from './modules/local/utils.nf' params(cleanupScript: cleanupScript,
																	 tmp: params.tmp, backupTmp: params.backupTmp,
																	 backupScratchHack: params.backupScratchHack)
include { } from './modules/local/step.0.generate_indexes.nf' params(
	refGenome: params.refGenome,
	refGenomeGtf: params.refGenomeGtf,
	erccFa: params.erccFa,
	erccGtf: params.erccGtf,
	starSjdbOverhang: params.starSjdbOverhang,
	publishDir: params.publishDir,
	hisatUseTranscript: params.hisatUseTranscript,
	tmp: params.tmp,
	backupTmp: params.backupTmp,
	backupScratchHack: params.backupScratchHack,
	nxfUnstageHack: params.nxfUnstageHack
)
include {
	download;
	check_direct_fastqs;
	filter_barcodes;
} from './modules/local/step.1.nf' params(
	parallelDownloads: params.parallelDownloads,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishFastqs,
	cleanupScript: cleanupScript,
	tmp: params.tmp,
	backupTmp: params.backupTmp,
	backupScratchHack: params.backupScratchHack,
	nxfUnstageHack: params.nxfUnstageHack
)
include {
	fastp;
} from './modules/local/step.2.nf' params(
	extraAdapters: params.extraAdapters,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishQCfiltered,
	cleanupScript: cleanupScript,
	tmp: params.tmp,
	backupTmp: params.backupTmp,
	backupScratchHack: params.backupScratchHack,
	nxfUnstageHack: params.nxfUnstageHack
)
include {
	star_counts;
	sort_bam;
	htseq_count;
} from './modules/local/step.readcounts.nf' params(
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishReadcounts,
	cleanupScript: cleanupScript,
	tmp: params.tmp,
	backupTmp: params.backupTmp,
	backupScratchHack: params.backupScratchHack,
	nxfUnstageHack: params.nxfUnstageHack
)
include {
	hisat2;
	ensure_hisat2_indexes;
} from './modules/local/step.3.hisat2.nf' params(
	hisatUseTranscript: params.hisatUseTranscript,
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishHisat,
	cleanupScript: cleanupScript,
	tmp: params.tmp,
	backupTmp: params.backupTmp,
	backupScratchHack: params.backupScratchHack,
	nxfUnstageHack: params.nxfUnstageHack
)
include {
	star;
	ensure_star_indexes;
} from './modules/local/step.4.star.nf' params(
	starUseSharedMem: params.starUseSharedMem,
	starThreadsSmall: params.starThreadsSmall,
	starThreadsLarge: params.starThreadsLarge,
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishStar,
	cleanupScript: cleanupScript,
	tmp: params.tmp,
	backupTmp: params.backupTmp,
	backupScratchHack: params.backupScratchHack,
	nxfUnstageHack: params.nxfUnstageHack
)
include {
	bowtie2;
	bowtie2_filter;
	ensure_bowtie2_indexes;
} from './modules/local/step.5.bowtie.nf' params(
	retainMixed: params.retainMixed, // retain mixed (xor) mate align cases
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishBowtie,
	cleanupScript: cleanupScript,
	tmp: params.tmp,
	backupTmp: params.backupTmp,
	backupScratchHack: params.backupScratchHack,
	nxfUnstageHack: params.nxfUnstageHack
)
include {
	dedup;
} from './modules/local/step.6.dedup.nf' params(
	percentLen: params.dedupPercentLen,
	minLen: params.dedupMinLen,
	pubishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishDedup,
	cleanupScript: cleanupScript,
	tmp: params.tmp,
	backupTmp: params.backupTmp,
	backupScratchHack: params.backupScratchHack,
	nxfUnstageHack: params.nxfUnstageHack
)
include {
	gsnap;
	gsnap_filter;
	gsnap_skip;
	ensure_gsnap_indexes;
} from './modules/local/step.7.gsnap.nf' params(
	retainMixed: params.retainMixed, // retain mixed (xor) mate align cases
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishGsnap,
	cleanupScript: cleanupScript,
	tmp: params.tmp,
	backupTmp: params.backupTmp,
	backupScratchHack: params.backupScratchHack,
	nxfUnstageHack: params.nxfUnstageHack
)
include {
	stats_csv;
} from './modules/local/stats.nf' params(
	skipHisat: params.skipHisat,
	tmp: params.tmp,
	backupTmp: params.backupTmp,
	backupScratchHack: params.backupScratchHack,
	nxfUnstageHack: params.nxfUnstageHack
)

CONTAINERS = ["docker", "singularity", "conda", "mamba"]
def container_usage() {
	return """Containerization options:
-profile docker           use docker containers to run commands when possible
-profile singularity      use singularity containers to run commands when possible
-profile conda            use conda packages to run commands when possible
-profile mamba            use conda packages, through mamba, to run commands when possible
"""
}

def join_by_id(ch1, ch2) {
	keyed_ch1_mates = ch1.map { meta, data -> [meta.id, meta, data] }
	keyed_ch2_mates = ch2.map { meta, data -> [meta.id, meta, data] }
	return keyed_ch1_mates.join(keyed_ch2_mates)
		.map { id, meta1, ch1_data, meta2, ch2_data -> [meta1, ch1_data, ch2_data] }
}

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
	// possibly check if container profiles are not active and we can't find a binary here
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
	// allow RunInfo csv (Run, size_MB)
	if (params.accessionList) {
		accessions = params.accessionList =~ /\.csv$/
			? Channel.fromPath(params.accessionList, type: 'file')
			.splitCsv( header: true )
			.map { row -> [row.Run.trim(), row.size_MB.toLong()] }
		: Channel.fromPath(params.accessionList, type: 'file')
			.splitText()
			.map { acc -> [acc.trim(), null] }
	} else {
		accessions = Channel.empty()
	}

	if (params.fastqPath) {
		try {
			direct_fastqs = Channel.fromFilePairs(params.fastqPath, size: -1, checkIfExists: true)
				.map { prefix, mates ->
					def size1 = mates[0].size()
					def new_meta = [id: prefix,
									sra_size: size1,
									size_MB: (mates[0].size() + mates[1]?.size()) / 1048576]
					[ new_meta, mates ]
				}
			direct_fastq_ids = direct_fastqs
				.map { meta, _ ->
					[ meta.id, null, true ]
				}

			// remove ids that exist in pre-dumped fastq path from accessions list
			// accessions = accessions
			// 	.join(direct_fastq_ids, remainder: true, by: 0)
			// 	.filter { key, s, v2 -> !v2 }
			// 	.map { key, MB, _ ->
			// 		[ [id: key, size_MB: MB, cleanup: "", cleanup_later: ""], key ]
			// 	}
			// 	.view()
		} catch (Exception e) {
			log.error "--fastq-path $params.fastqPath is not folders of fastq?\n$e"
			exit(1)
		}
	} else {
		direct_fastqs = Channel.empty()
		accessions = accessions
			.map { key, MB ->
				[ [id: key, size_MB: MB, cleanup: "", cleanup_later: ""], key ]
			}
	}
	
	// download SRAs by remaining accessions
	download(accessions).mates
		.map { meta, fastq, reads, sra_size, median1, median2, count, fsize ->
			def new_meta = [id: meta.id,
							reads: count.toLong(),
							sra_size: sra_size.toLong(),
							readlen: median1.toLong(),
							readlen_2: median2 != "" ? median2.toLong() : "",
							fastq_size: fsize.toLong(),
							single_end: fastq.size() != 2,
							size_MB: meta.size_MB,
							cleanup: "",
							cleanup_later: "${fastq.join(' ')}"]
			[ new_meta, fastq ]
		}
		.view()
		.set { download_result }

	direct_fastqs.view()
	n_direct_fastqs = check_direct_fastqs(direct_fastqs) // , merging any existing fastqs
		.map {
			meta, fastq, lines ->
			def new_meta = [id: meta.id,
							reads: lines.toLong() / 4,
							sra_size: meta.sra_size.toLong(),
							size_MB: meta.size_MB,
							cleanup: "",
							cleanup_later: "${fastq.join(' ')}"]
			[ new_meta, fastq ]
		}
		.view()

	// heuristic filter scRNAseq barcode files and add metadata for direct path
	filter_barcodes(n_direct_fastqs)
		.map { meta, fastq, median1, median2, count, fsize ->
			def new_meta = meta.clone()
			new_meta.reads = count.toLong()
			new_meta.readlen = median1.toLong()
			new_meta.readlen_2 = median2 != "" ? median2.toLong() : ""
			new_meta.fastq_size = fsize.toLong()
			new_meta.single_end = fastq.size() != 2
			new_meta.size_MB = meta?.size_MB
			new_meta.cleanup = meta.cleanup_later
			new_meta.cleanup_later = "${fastq.join(' ')}"
			[ new_meta, fastq ]
		}
		.mix(download_result)
		.branch { // empty/insignificant runs (by trimming, qc, or host mapping) should drop out
			ok: { meta, fastq -> 
				meta.single_end ? file(fastq[0]).size() > 132 : (fastq.size() == 2) &&
					file(fastq[0]).size() > 132 && file(fastq[1]).size() > 132
			}(it)
			dropouts: true
		}
		.set { filter_barcodes_result }


	// step 2: adapter trimming & filtering
	fastp(filter_barcodes_result.ok)
	fastp.out.mates
		.map { meta, fastq, fastp_reads_after -> {
				def new_meta = meta.clone()
				new_meta.cleanup = meta.cleanup_later
				if (!params.skipHostCounts) {
					new_meta.cleanup_later = "" // fastp is before branch so needs special cleanup
					new_meta.qc_cleanup = "${fastq.join(' ')}"
				} else {
					new_meta.cleanup_later = "${fastq.join(' ')}"
					new_meta.qc_cleanup = ""
				}
				new_meta.fastp_reads_after = fastp_reads_after.toLong()
				[ new_meta, fastq ]
			}
		}
		.branch { // empty/insignificant runs (by trimming, qc, or host mapping) should drop out)
			ok: { meta, fastq -> 
				meta.single_end ? file(fastq[0]).size() > 132 : (fastq.size() == 2) &&
					file(fastq[0]).size() > 132 && file(fastq[1]).size() > 132
			}(it)
			
			dropouts: true
		}
		.set { fastp_result }


	// host read counts path
	if (!params.skipHostCounts) {
		star_counts(fastp_result.ok, star_indexes2).bam
			.map { meta, bam ->
				def new_meta = meta.clone()
				// cleanup didn't happen in this step
				new_meta.cleanup = "${meta.cleanup} ${meta.cleanup_later}"
				new_meta.cleanup_later = "${bam.toString()}"
				[ new_meta, bam ]
			}
			.set { starcounts_result }
		sort_bam(starcounts_result).bam
			.map { meta, bam ->
				def new_meta = meta.clone()
				// cleanup didn't happen in this step
				new_meta.cleanup = "${meta.cleanup} ${meta.cleanup_later}"
				new_meta.cleanup_later = "${bam.toString()}"
				[ new_meta, bam ]
			}
			.set { sortbam_result }

		htseq_count(sortbam_result, file(params.refGenomeGtf))
			.set { count_result }
	}

	// step 3: hisat2
	if (!params.skipHisat) {
		hisat2(fastp_result.ok, hisat2_indexes)
		hisat2.out.mates
			.map { meta, mates ->
				m = meta.clone(); m.cleanup = m.cleanup_later; m.cleanup_later = "${mates.join(' ')}"
				[ m, mates ]
			}
			.set { unmapped_reads_1 }
	} else {
		unmapped_reads_1 = fastp_result.ok
	}

	// step 4: STAR
	star(unmapped_reads_1, star_indexes).mates
		.map { meta, mates ->
			m = meta.clone(); m.cleanup = m.cleanup_later; m.cleanup_later = "${mates.join(' ')}"
			[ m, mates ]
		}
		.set { star_result }

	// clean up the pre-branch checkpoint
	if (!params.skipHostCounts) {
		cleanup_branched(join_by_id(star_result, sortbam_result))
	}
	
	// step 5: bowtie2
	bowtie2(star_result, bowtie2_indexes).sam
		.map { meta, sam ->
			m = meta.clone(); m.cleanup = m.cleanup_later; m.cleanup_later = "${sam.toString()}"
			[ m, sam ]
		}
		.set { bowtie2_result }
	bowtie2_filter_input = join_by_id(bowtie2_result, star.out.mates)
	bowtie2_filter(bowtie2_filter_input)

	// step 6: deduplication
	bowtie2_filter.out.mates
		.map { meta, mates ->
			m = meta.clone(); m.cleanup = m.cleanup_later; m.cleanup_later = "${mates.join(' ')}"
			[ m, mates ]
		}
		.branch { // empty/insignificant runs (by trimming, qc, or host mapping) should drop out)
			ok: { meta, fastq ->
				meta.single_end ? file(fastq[0]).size() > 132 : (fastq.size() == 2) &&
					file(fastq[0]).size() > 132 && file(fastq[1]).size() > 132
			}(it)
			dropouts: true
		}
		.set { bowtie2_filtered_result }

	join_by_id(bowtie2_filtered_result.ok, star.out.stats)
		.set { dedup_input }
	dedup(dedup_input)
	// don't add dedup to cleanup list

	// step 7: gsnap
	gsnap(dedup.out.mates, gsnap_indexes)
	gsnap.out.sam
		.map { meta, sam ->
			m = meta.clone(); m.cleanup = m.cleanup_later; m.cleanup_later = "${sam.toString()}"
			[ m, sam ]
		}
		.branch {
			ok: { meta, sam -> file(sam).size() > 0 }(it)
			fails: true
		}
		.set { gsnap_result }
	gsnap_filter_input = join_by_id(gsnap_result.ok, dedup.out.mates)
	gsnap_filter(gsnap_filter_input)

	// if gsnap fails a run for any non-oom reason, keep the nonhost reads from before that.
	dedup.out.mates
		.join(gsnap.out.sam, by: [0], remainder: true)
		.filter { meta, dedup_mates, gsnap -> !gsnap || file(gsnap).size() == 0 }
		.map { meta, dedup_mates, gsnap_null -> [ meta, dedup_mates ] }
		.set { gsnap_skips }
	gsnap_skip(gsnap_skips)

	// rekey in preparation for join
	// questionable to use prefetch / fastq_dump stats b/c not present for direct fastq
	// but we DO
	//sra_stats = sra.out.stats.map { meta, sra -> [ meta.id, sra ] }
	//fastq_stats = fastqs.out.stats.map { meta, fastq -> [ meta.id, sra ] }
	meta_stats    = filter_barcodes_result.ok.map { meta, fastq -> [ meta.id, meta ] }
	fastp_stats   = fastp.out.stats_txt.map       { meta, stats -> [ meta.id, stats ] }
	if (!params.skipHisat) {
		hisat2_stats  = hisat2.out.stats.map  { meta, stats -> [ meta.id, stats ] }
	} else {
		hisat2_stats = meta_stats.map { id, _ -> [ id, "na" ] }
	}
	star_stats    = star.out.stats.map    { meta, stats -> [ meta.id, stats ] }
	bowtie2_stats = bowtie2_filter.out.stats.map { meta, stats -> [ meta.id, stats ] }
	dedup_stats   = dedup.out.stats.map   { meta, stats -> [ meta.id, stats ] }
	gsnap_stats   = gsnap_filter.out.stats.map { meta, stats -> [ meta.id, stats, "yes" ] }
		.concat (gsnap_skip.out.mates.map { meta, mates -> [ meta.id, null, "no" ] })

	all_stats = meta_stats.join(fastp_stats)
		.join(hisat2_stats)
		.join(star_stats)
		.join(bowtie2_stats)
		.join(dedup_stats)
		.join(gsnap_stats)

	stats_csv(all_stats)
		.collectFile(name: "stats-${params.timestamp}.csv", keepHeader: true, skip: 1, storeDir: "${params.publishDir}/stats/")
}
