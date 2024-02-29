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

include {
	bowtie2_generate_indexes;
} from './step.0.generate_indexes.nf' params(
	publishDir: params.publishDir,
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
	tuple val(idx_basename), path("bowtie2_index/*")

	output:
	tuple val(meta), path("bowtie2.sam"), emit: sam

	script:
	// can we safely enable --no-discordant here?
	//def index_base = file(index_dir).listFiles()[0].getSimpleName()
	def BOWTIE2_CMD = """bowtie2 --quiet --very-sensitive-local -p $task.cpus \
		--rg-id na --rg LB:na --rg SM:na --rg PL:na --rg PU:na \
		-x bowtie2_index/${idx_basename} """
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
	tuple val(meta), path("bowtie2.sam")

	output:
	tuple val(meta), path("bowtie2.unmapped.names.txt"), emit: names
	tuple val(meta), path("bowtie2.stats.txt"), emit: stats

	// +discordant would be (flag.paired && !flag.proper_pair)
	script:
	def cond = params.retainMixed ?
		'flag.unmap || (flag.paired && flag.munmap)' :
		'(!flag.paired && flag.unmap) || (flag.paired && flag.unmap && flag.munmap)'
	//def FILTER_CMD = "LC_ALL=C grep -A3 --color=never -wFf bowtie2.unmapped.names.txt | sed '/^--\$/d'"
	"""
	samtools view -@ ${task.cpus} bowtie2.sam | cut -f2 | sort | uniq -c > bowtie2.stats.txt
	samtools view -@ ${task.cpus} -e '$cond' bowtie2.sam | cut -f1 > bowtie2.unmapped.names.txt
	"""
}

process bowtie2_filter_by_names {
	label "fastq_filter"

	input:
	tuple val(meta), path(names), path(mategz)

	def SUFFIX = "filteredbyBT.fastq"
	output:
	tuple val(meta), path("Unmapped.out.mate?.${SUFFIX}.gz"), emit: mates

	script:
	def FILTER_CMD = "LC_ALL=C fastq-namefilter $names -"
	if (!meta.single_end)
	"""
	gunzip -c ${mategz[0]} | $FILTER_CMD > Unmapped.out.mate1.${SUFFIX}
	gunzip -c ${mategz[1]} | $FILTER_CMD > Unmapped.out.mate2.${SUFFIX}
	gzip -nc Unmapped.out.mate1.${SUFFIX} > Unmapped.out.mate1.${SUFFIX}.gz.staging
	gzip -nc Unmapped.out.mate2.${SUFFIX} > Unmapped.out.mate2.${SUFFIX}.gz.staging
	mv Unmapped.out.mate1.${SUFFIX}.gz.staging Unmapped.out.mate1.${SUFFIX}.gz
	mv Unmapped.out.mate2.${SUFFIX}.gz.staging Unmapped.out.mate2.${SUFFIX}.gz
	"""
	else if (meta.single_end)
	"""
	gunzip -c ${mategz} | $FILTER_CMD > Unmapped.out.mate1.${SUFFIX}
	gzip -nc Unmapped.out.mate1.${SUFFIX} > Unmapped.out.mate1.${SUFFIX}.gz.staging
	mv Unmapped.out.mate1.${SUFFIX}.gz.staging Unmapped.out.mate1.${SUFFIX}.gz
	"""
}


def ensure_bowtie2_indexes(ref_indexes,
						   ref_genome, ercc) {
	if (ref_indexes
		&& (indexes = file(ref_indexes))
		&& indexes.exists()) {
		return [indexes.listFiles()[0].getSimpleName(), file(indexes.resolve("*.{bt2,bt2l}"))]
	} else {
		indexes = bowtie2_generate_indexes(file(ref_genome),
										   file(ercc))
	}
	return indexes
}