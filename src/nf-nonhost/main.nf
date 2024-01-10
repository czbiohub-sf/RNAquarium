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
params.hisatUseTranscript = false
params.help = false

params.genomeSize = null // must be filled
params.starRefIndexesErcc = null // "Danio_rerio.GRCz11.108.ERCC"
params.starRefIndexes = null // "Danio_rerio.GRCz11.108"
params.hisatRefIndexes = null
params.bowtieRefIndexes = null
params.gmapRefIndexes = null
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
// consider removing these options, which are potentially dangerous for
// public-exposed pipelines and have questionable utility even for y!!
// params.starIndexGenOptions = ""
// params.hisatIndexGenOptions = ""
// params.bowtieIndexGenOptions = ""
// params.gsnapIndexGenOptions = ""
// params.sraPrefetchOptions = ""
// params.fastqDumpOptions = ""
// params.fastpOptions = ""
// params.priceOptions = ""
// params.starCountOptions = "" // host read counts pipeline
// params.samtoolsSortOptions = ""
// params.htseqCountOptions = ""
// params.hisatOptions = "" // nonhost pipeline
// params.starOptions = ""
// params.bowtieOptions = ""
// params.gsnapOptions = ""

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'
include { } from './modules/step.0.generate_indexes.nf' params(
	publishDir: params.publishDir,
	hisatUseTranscript: params.hisatUseTranscript
)
include {
	prefetch;
	fastq_dump;
	check_direct_fastqs;
	filter_barcodes;
} from './modules/step.1.nf' params(
	parallelDownloads: params.parallelDownloads,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishFastqs,
)
include {
	fastp;
	priceseqfilter;
} from './modules/step.2.nf' params(
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishPricefiltered,
)
include {
	star_counts;
	sort_bam;
	htseq_count;
} from './modules/step.readcounts.nf' params(
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishReadcounts,
)
include {
	hisat2;
	ensure_hisat2_indexes;
} from './modules/step.3.hisat2.nf' params(
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishHisat,
)
include {
	star;
	ensure_star_indexes;
} from './modules/step.4.star.nf' params(
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
} from './modules/step.5.bowtie.nf' params(
	retainMixed: params.retainMixed, // retain mixed (xor) mate align cases
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishBowtie,
)
include {
	dedup;
} from './modules/step.6.dedup.nf' params(
	percentLen: params.dedupPercentLen,
	minLen: params.dedupMinLen,
	pubishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishDedup,
)
include {
	gsnap;
	process_gsnap_sam;
	gsnap_filter_by_names;
	ensure_gsnap_indexes;
} from './modules/step.7.gsnap.nf' params(
	retainMixed: params.retainMixed, // retain mixed (xor) mate align cases
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	publishIntermediate: params.publishIntermediate && params.publishGsnap,
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
	gsnap_indexes = ensure_gsnap_indexes(params.gmapRefIndexes,
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
	fastp(fastqs_2)
	filtered_fastqs = priceseqfilter(fastp.out.mates).mates

	// host read counts path
	if (!params.skipHostCounts) {
		star_counts(filtered_fastqs, star_indexes2)
		sort_bam(star_counts.out.bam)
		count = htseq_count(sort_bam.out.bam, file(params.refGenomeGtf))
	}

	// step 3: hisat2
	if (!params.skipHisat) {
		hisat2(filtered_fastqs, hisat2_indexes)
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
	dedup(bowtie2_filter_by_names.out.mates.join(star.out.stats))

	// step 7: gsnap
	gsnap(dedup.out.mates, gsnap_indexes)
	process_gsnap_sam(gsnap.out.sam)
	gsnap_filter_by_names(process_gsnap_sam.out.names.join(dedup.out.mates))

	// rekey in preparation for join
	// questionable to use prefetch / fastq_dump stats bc not present for direct fastq
	//sra_stats = sra.out.stats.map { meta, sra -> [ meta.id, sra ] }
	//fastq_stats = fastqs.out.stats.map { meta, fastq -> [ meta.id, sra ] }
	meta_stats    = fastqs_2.map          { meta, fastq -> [ meta.id, meta ] }
	fastp_stats   = fastp.out.stats.map   { meta, stats -> [ meta.id, stats ] }
	price_stats   = priceseqfilter.out.stats.map { meta, stats -> [ meta.id, stats ] }
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
		.collectFile(name: "stats-${params.timestamp}.csv", keepHeader: true, skip: 1, storeDir: "$PWD")
		.view()
}

process stats_csv {
	cache = false

	input:
	tuple val(idx), val(meta), file("fastp_stats.txt"), val(price_stats), file("hisat2_stats.txt"), file("star_stats.txt"), file("bowtie2_stats.txt"), val(dedup_stats), file("gsnap_stats.txt")

	output:
	stdout

	script:
	def HISAT_UNPAIRED = '''hisat_before=\$(sed -n '/reads; of these:/{;p;}' hisat2_stats.txt | cut -f1 -d' ')
	hisat_unaligned=\$(sed -n '/aligned 0 times/{;p;}' hisat2_stats.txt | cut -f5 -d' ')
	hisat_aligned_unique=\$(sed -n '/aligned exactly 1 time/{;p;}' hisat2_stats.txt | cut -f5 -d' ')
	hisat_multialign=\$(sed -n '/aligned >1 times/{;p;}' hisat2_stats.txt | cut -f5 -d' ')
	hisat_discordant=''
	printf "%s,%s,%s,%s,%s," "\$hisat_before" "\$hisat_unaligned" "\$hisat_aligned_unique" "\$hisat_multialign" "\$hisat_discordant"
	'''
	def HISAT_PAIRED = '''hisat_before=\$(sed -n '/reads; of these:/{;p;}' hisat2_stats.txt | cut -f1 -d' ')
	hisat_unaligned=\$(sed -n '/pairs aligned concordantly 0 times/{;p;}' hisat2_stats.txt | cut -f5 -d' ')
	hisat_aligned_unique=\$(sed -n '/aligned concordantly exactly 1 time/{;p;}' hisat2_stats.txt | cut -f5 -d' ')
	hisat_multialign=\$(sed -n '/aligned concordantly >1 times/{;p;}' hisat2_stats.txt | cut -f5 -d' ')
	hisat_discordant=\$(sed -n '/aligned discordantly 1 time/{;p;}' hisat2_stats.txt | cut -f5 -d' ')
	printf "%s,%s,%s,%s,%s," "\$hisat_before" "\$hisat_unaligned" "\$hisat_aligned_unique" "\$hisat_multialign" "\$hisat_discordant"
	'''
	def PARSE_HISAT = params.skipHisat ? ",,,,," : (meta.single_end ? HISAT_UNPAIRED : HISAT_PAIRED)
	"""
	#       id single_end     reads       readlen
	printf "id,single_end,starting_reads,median_len,"
	#
	printf "fastp_reads_before,fastp_reads_after,fastp_reads_too_short,fastp_reads_trimmed,"
	#                               %1 - %2
	printf "price_reads_before,price_reads_after,"
	printf "hisat2_reads_before,hisat2_unaligned,hisat2_aligned_unique,hisat2_multialign,hisat2_discordant,"
	printf "star_reads_before,star_avg_len,star_aligned_unique,star_multialign,star_unaligned,star_too_short,"
	printf "bowtie2_reads_before,bowtie2_aligned,bowtie2_multialign,bowtie2_aligned_unique,bowtie2_unaligned,bowtie2_mixed,"
	printf "dedup_reads_before,dedup_reads_after,"
	printf "gsnap_reads_before,gsnap_aligned,gsnap_multialign,gsnap_aligned_unique,gsnap_unaligned,gsnap_mixed\n"

	printf "${idx},${meta.single_end},${meta.reads},${meta.readlen},"

	# FASTP
	fastp_before=\$(sed -n '/Read1 before filtering:/{;n;p;}' fastp_stats.txt | cut -f3 -d' ')
	fastp_after=\$(sed -n '/Read1 after filtering:/{;n;p;}' fastp_stats.txt | cut -f3 -d' ')
	fastp_short=\$(sed -n '/reads failed due to too short:/{;p;}' fastp_stats.txt | cut -f7 -d' ') 
	fastp_trimmed=\$(sed -n '/reads with adapter trimmed:/{;p;}' fastp_stats.txt | cut -f5 -d' ')
	printf "%s,%s,%s,%s," \$fastp_before \$fastp_after \$fastp_short \$fastp_trimmed

	# PRICE
	# price outputs \b to non-interactive output
	price_total=\$(echo "${price_stats.split("\n")[3].split("\b")[-1].split("/")[0]}" | grep -o "[0-9]\\+")
	price_removed=\$(echo "${price_stats.split("\n")[3].split("\b")[-1].split("/")[1]}" | grep -o "[0-9]\\+")
	price_after=\$(bc <<< "\$price_total - \$price_removed")
	printf "%d,%d," "\$price_total" "\$price_after"

	# HISAT2
	${PARSE_HISAT}

	# STAR
	star_before=\$(grep "Number of input reads |" star_stats.txt | grep -o "[0-9]\\+")
	star_avg=\$(grep "Average input read length |" star_stats.txt | grep -o "[0-9]\\+")
	star_unique=\$(grep "Uniquely mapped reads number |" star_stats.txt | grep -o "[0-9]\\+")
	star_multialign=\$(grep "Number of reads mapped to multiple loci |" star_stats.txt | grep -o "[0-9]\\+")
	star_unaligned=\$(bc <<< "\$star_before - \$star_unique - \$star_multialign")
	star_short=\$(grep "Number of reads unmapped: too short |" star_stats.txt | grep -o "[0-9]\\+")
	printf "%s,%s,%s,%s,%s,%s," \$star_before \$star_avg \$star_unique \$star_multialign \$star_unaligned \$star_short

	# BOWTIE2
	bowtie2_total=0
	bowtie2_multi=0
	bowtie2_aligned=0
	bowtie2_unaligned=0
	bowtie2_mixed=0
	while IFS=' ' read -r count sambits
	do
		# ignore supplementary alignments and second-in-pair
		if [[ \$(( \$sambits & 0x80 || \$sambits & 0x800 )) -ne 0 ]]; then
			continue
		fi
		bowtie2_total=\$(( \$bowtie2_total + \$count ))
		# track secondary mappings but don't add to other counts
		if [[ \$(( "\$sambits" & 0x100 )) -ne 0 ]]; then
			bowtie2_multi=\$(( "\$bowtie2_multi" + \$count ))
			continue
		fi
		if [[ \$(( (\$sambits & 0x5)==0x4 || (\$sambits & 0xD)>0x4 )) -ne 0 ]]; then
			bowtie2_unaligned=\$(( \$bowtie2_unaligned + \$count ))
		fi
		if [[ \$(( ! (\$sambits & 0x4 || \$sambits & 0x8) )) -ne 0 ]]; then
			bowtie2_aligned=\$(( \$bowtie2_aligned + \$count ))
		fi
		if [[ \$(( (\$sambits & 0xD)==0x5 || (\$sambits & 0xD)==0x9  )) -ne 0 ]]; then
			bowtie2_mixed=\$(( \$bowtie2_mixed + \$count ))
		fi
	done < bowtie2_stats.txt
	bowtie2_unique=\$(( \$bowtie2_aligned - \$bowtie2_multi ))
	printf "%s,%s,%s,%s,%s,%s," \$bowtie2_total \$bowtie2_aligned \$bowtie2_multi \$bowtie2_unique \$bowtie2_unaligned \$bowtie2_mixed

	# DEDUP
	dedup_before=\$(echo "$dedup_stats" | grep "total reads:"  | grep -o "[0-9]\\+")
	dedup_after=\$(echo "$dedup_stats" | grep "unique reads:"  | grep -o "[0-9]\\+")
	printf "%s,%s," "\$dedup_before" "\$dedup_after"

	# GSNAP
	gsnap_total=0
	gsnap_multi=0
	gsnap_aligned=0
	gsnap_unaligned=0
	gsnap_mixed=0
	while IFS=' ' read -r count sambits
	do
		# ignore supplementary alignments and second-in-pair
		if [[ \$(( \$sambits & 0x80 || \$sambits & 0x800 )) -ne 0 ]]; then
			continue
		fi
		gsnap_total=\$(( \$gsnap_total + \$count ))
		# track secondary mappings but don't add to other counts
		if [[ \$(( "\$sambits" & 0x100 )) -ne 0 ]]; then
			gsnap_multi=\$(( "\$gsnap_multi" + \$count ))
			continue
		fi
		if [[ \$(( (\$sambits & 0x5)==0x4 || (\$sambits & 0xD)>0x4 )) -ne 0 ]]; then
			gsnap_unaligned=\$(( \$gsnap_unaligned + \$count ))
		fi
		if [[ \$(( ! (\$sambits & 0x4 || \$sambits & 0x8) )) -ne 0 ]]; then
			gsnap_aligned=\$(( \$gsnap_aligned + \$count ))
		fi
		if [[ \$(( (\$sambits & 0xD)==0x5 || (\$sambits & 0xD)==0x9  )) -ne 0 ]]; then
			gsnap_mixed=\$(( \$gsnap_mixed + \$count ))
		fi
	done < gsnap_stats.txt
	gsnap_unique=\$(( \$gsnap_aligned - \$gsnap_multi ))
	printf "%s,%s,%s,%s,%s,%s\n" \$gsnap_total \$gsnap_aligned \$gsnap_multi \$gsnap_unique \$gsnap_unaligned \$gsnap_mixed
	"""
}