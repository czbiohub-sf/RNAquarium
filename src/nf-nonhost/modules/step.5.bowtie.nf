#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.publishDir = "$PWD"
params.publishIntermediate = true
params.genomeSize = null
params.bowtieIndex = null

params.refGenome = "Danio_rerio.GRCz11.dna_sm.primary_assembly.fa"
params.refGenomeGtf = "Danio_rerio.GRCz11.108.gtf"
params.erccFa = "ERCC92.fa"
params.erccGtf = "ERCC92.gtf"

params.retainMixed = true
params.bowtieOptions = ""
params.bowtieIndexGenOptions = ""

include {
	bowtie2_generate_indexes;
} from './step.0.generate_indexes.nf' params(
	publishDir: params.publishDir,
	bowtieIndexGenOptions: params.bowtieIndexGenOptions
)

params.metaIn = 'step_4_sheet.csv'
params.metaOut = 'step_5_sheet.csv'
include {
	LOAD_METASHEET;
	SAVE_METASHEET;
} from './utils.nf'

 process bowtie2 {
	label 'bowtie2'

	input:
	tuple val(meta), path(mategz)
	path index_dir

	output:
	tuple val(meta), path("bowtie2.sam"), emit: sam

	script:
	// can we safely enable --no-discordant here?
	def index_base = file(index_dir).listFiles()[0].getBaseName(2)
	def BOWTIE2_CMD = """bowtie2 --quiet --very-sensitive-local -p $task.cpus \
		--rg-id na --rg LB:na --rg SM:na --rg PL:na --rg PU:na \
		-x $index_dir/$index_base $params.bowtieOptions"""
	if (!meta.single_end)
	"""
	${BOWTIE2_CMD} -1 ${mategz[0]} -2 ${mategz[1]} -S bowtie2.sam.staging
	mv bowtie2.sam.staging bowtie2.sam
	"""
	else if (meta.single_end)
	"""
	${BOWTIE2_CMD} -U ${mategz} -S bowtie2.sam.staging
	mv bowtie2.sam.staging bowtie2.sam
	"""
}

process process_bowtie2_sam {
	label 'samtools'

	input:
	tuple val(_), path("bowtie2.sam")
	tuple val(meta), path(mategz)

	output:
	tuple val(meta), path("*.mate?.filteredbyBT.gz")
	tuple val(meta), path("bowtie2.stats.txt"), emit: stats

	// +discordant would be (flag.paired && !flag.proper_pair)
	script:
	def cond = params.retainMixed ?
		'flag.unmap || (flag.paired && flag.munmap)' :
		'(!flag.paired && flag.unmap) || (flag.paired && flag.unmap && flag.munmap)'
	def PREFILTER = """samtools view -b bowtie2.sam > bowtie2.bam
	samtools view bowtie2.bam | cut -f2 | sort | uniq -c > bowtie2.stats.txt
	samtools view -e '$cond' bowtie2.bam | cut -f1 | sed 's/^/@/' > bowtie2.unmapped.names.txt"""
	def FILTER_CMD = "grep -A3 --color=never -wFf bowtie2.unmapped.names.txt | sed '/^--\$/d'"
	if (!meta.single_end)
	"""
	echo $meta
	${PREFILTER}
	gunzip -c ${mategz[0]} | $FILTER_CMD > Unmapped.out.mate1.filteredbyBT
	gunzip -c ${mategz[1]} | $FILTER_CMD > Unmapped.out.mate2.filteredbyBT
	gzip -c Unmapped.out.mate1.filteredbyBT > Unmapped.out.mate1.filteredbyBT.gz.staging
	gzip -c Unmapped.out.mate2.filteredbyBT > Unmapped.out.mate2.filteredbyBT.gz.staging
	mv Unmapped.out.mate1.filteredbyBT.gz.staging Unmapped.out.mate1.filteredbyBT.gz
	mv Unmapped.out.mate2.filteredbyBT.gz.staging Unmapped.out.mate2.filteredbyBT.gz
	"""
	else if (meta.single_end)
	"""
	${PREFILTER}
	gunzip -c ${mategz} | $FILTER_CMD > Unmapped.out.mate1.filteredbyBT
	gzip -c Unmapped.out.mate1.filteredbyBT > Unmapped.out.mate1.filteredbyBT.gz.staging
	mv Unmapped.out.mate1.filteredbyBT.gz.staging Unmapped.out.mate1.filteredbyBT.gz
	"""
}


def ensure_bowtie2_indexes(ref_indexes,
						   ref_genome, ercc) {
	if (ref_indexes) {
		indexes = file(ref_indexes)
	} else {
		indexes = bowtie2_generate_indexes(file(ref_genome),
										   file(ercc))
	}
	return indexes
}