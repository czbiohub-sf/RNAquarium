#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.publishDir = "$PWD"
params.publishIntermediate = true
params.genomeSize = null
params.maxMismatch = 0.3

params.metaIn = 'step_6_sheet.csv'
params.metaOut = 'step_7_sheet.csv'
include {
	LOAD_METASHEET;
	SAVE_METASHEET;
} from './utils.nf'

include {
	gsnap_generate_indexes;
} from './step.0.generate_indexes.nf' params(
	genomeSize: params.genomeSize,
	publishDir: params.publishDir,
	bowtieIndexGenOptions: params.bowtieIndexGenOptions
)

process gsnap {
	label 'gmap'

	input:
	tuple val(meta), path(mategz)
	path index_dir

	output:
	tuple val(meta), path("gsnap_out.sam"), emit: sam

	script:
	def index_name = file(index_dir).getName()
	def gsnap_gmap_bin = params.genomeSize < 2**32 ? "gsnap" : "gmap"
	def GSNAP_CMD = """${gsnap_gmap_bin} -A sam \
		-N 1 `# find novel splice sites` \
		--batch=2 \
		--use-shared-memory=0 \
		--npaths=1 `# maximum paths to print` \
		--ordered -t ${task.cpus} \
		--max-mismatches $params.maxMismatch \
		-D $index_dir -d $index_name \
		-o gsnap_out.sam.staging """
	if (!meta.single_end)
	"""
	gunzip -kcd ${mategz[0]} > mate1.fastq
	gunzip -kcd ${mategz[1]} > mate2.fastq
	${GSNAP_CMD} mate1.fastq mate2.fastq
	mv gsnap_out.sam.staging gsnap_out.sam
	"""
	else if (meta.single_end)
	"""
	gunzip -kcd ${mategz} > mate1.fastq
	${GSNAP_CMD} mate1.fastq
	mv gsnap_out.sam.staging gsnap_out.sam
	"""
}

process process_gsnap_sam {
	label 'samtools'

	input:
	tuple val(meta), path(gsnap_sam)

	output:
	tuple val(meta), path("gsnap.unmapped.names.txt"), emit: names
	tuple val(meta), path("gsnap.stats.txt"), emit: stats

	script:
	def cond = params.retainMixed ?
		'flag.unmap || (flag.paired && flag.munmap)' :
		'(!flag.paired && flag.unmap) || (flag.paired && flag.unmap && flag.munmap)'
	//def FILTER_CMD = "LC=ALL grep -A3 --color=never -wFf gsnap.unmapped.names.txt | sed '/^--\$/d'"
	"""
	samtools view -@ ${task.cpus} $gsnap_sam | cut -f2 | sort | uniq -c > gsnap.stats.txt
	samtools view -@ ${task.cpus} -e '$cond' $gsnap_sam | cut -f1  > gsnap.unmapped.names.txt
	"""
}

process gsnap_filter_by_names {
	label "fastq_filter"

	input:
	tuple val(meta), path(names), path(mategz)

	def SUFFIX = "filteredbyBT.dedup.gsnapFiltered.fastq"
	output:
	tuple val(meta), path("Unmapped.out.mate?.${SUFFIX}.gz"), emit: mates

	script:
	def FILTER_CMD = "LC_ALL=C fastq-namefilter $names -"
	if (!meta.single_end)
	"""
	gunzip -c ${mategz[0]} | $FILTER_CMD > Unmapped.out.mate1.${SUFFIX}
	gunzip -c ${mategz[1]} | $FILTER_CMD > Unmapped.out.mate2.${SUFFIX}
	gzip -c Unmapped.out.mate1.${SUFFIX} > Unmapped.out.mate1.${SUFFIX}.gz.staging
	gzip -c Unmapped.out.mate2.${SUFFIX} > Unmapped.out.mate2.${SUFFIX}.gz.staging
	mv Unmapped.out.mate1.${SUFFIX}.gz.staging Unmapped.out.mate1.${SUFFIX}.gz
	mv Unmapped.out.mate2.${SUFFIX}.gz.staging Unmapped.out.mate2.${SUFFIX}.gz
	"""
	else if (meta.single_end)
	"""
	gunzip -c ${mategz} | $FILTER_CMD > Unmapped.out.mate1.${SUFFIX}
	gzip -c Unmapped.out.mate1.${SUFFIX} > Unmapped.out.mate1.${SUFFIX}.gz.staging
	mv Unmapped.out.mate1.${SUFFIX}.gz.staging Unmapped.out.mate1.${SUFFIX}.gz
	"""
}

def ensure_gsnap_indexes(ref_indexes,
						 ref_genome, ercc) {
	if (ref_indexes
		&& (indexes = file(ref_indexes))
		&& indexes.exists()) {
	} else {
		indexes = gsnap_generate_indexes(file(ref_genome),
										 file(ercc))
	}
	return indexes
}