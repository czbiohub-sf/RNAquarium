#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.publishDir = "$PWD"
params.publishIntermediate = true
params.genomeSize = null
params.refIndexes = null
params.cleanupScript = ""

params.refGenome = "Danio_rerio.GRCz11.dna_sm.primary_assembly.fa"
params.refGenomeGtf = "Danio_rerio.GRCz11.108.gtf"
params.erccFa = "ERCC92.fa"
params.erccGtf = "ERCC92.gtf"

include {
	hisat2_generate_indexes;
} from './step.0.generate_indexes.nf' params(
	publishDir: params.publishDir,
)

params.metaIn = 'step_2_sheet.csv'
params.metaOut = 'step_3_sheet.csv'
include {
	LOAD_METASHEET;
	SAVE_METASHEET;
} from './utils.nf'

process hisat2 {
	label 'hisat2'
	publishDir "$params.publishDir/hisat2_out/$meta.id", enabled: params.publishIntermediate

	input:
	tuple val(meta), path(fqgz, arity: '1..2')
	tuple val(idx_basename), path("hisat2_index/*")

	output:
	tuple val(meta), path("?E/Unmapped.out.mate?.gz", arity: '1..2'), emit: mates
	tuple val(meta), path("stats.txt"), emit: stats
	tuple val(meta), path("metrics.txt"), emit: hisat2_debug
	
	script:
	def HISAT2_CMD = """hisat2 --met-file metrics.txt --summary-file stats.txt -p $task.cpus -k 1 -S /dev/null \
		-x hisat2_index/${idx_basename} """
	if (!meta.single_end) """
	mkdir -p PE
	${HISAT2_CMD} \
		-1 ${fqgz[0]} -2 ${fqgz[1]} \
		--un-conc-gz PE/Unmapped.out.mate%.gz.staging
	mv PE/Unmapped.out.mate1.gz.staging PE/Unmapped.out.mate1.gz
	mv PE/Unmapped.out.mate2.gz.staging PE/Unmapped.out.mate2.gz

	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
	else if (meta.single_end) """
	mkdir -p SE
	${HISAT2_CMD} \
		-U ${fqgz} \
		--un-gz SE/Unmapped.out.mate1.gz.staging
	mv SE/Unmapped.out.mate1.gz.staging SE/Unmapped.out.mate1.gz

	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
}

def ensure_hisat2_indexes(ref_indexes,
						  ref_genome, ref_genome_gtf, ercc, ercc_gtf) {
	if (ref_indexes
		&& (indexes = file(ref_indexes))
		&& indexes.exists()) {
		return [indexes.listFiles()[0].getSimpleName(), file(indexes.resolve("*.{ht2,ht2l}")) ]
	} else {
		indexes = hisat2_generate_indexes(file(ref_genome),
										  file(ref_genome_gtf),
										  file(ercc),
										  file(ercc_gtf))
	}
	return indexes
}


workflow {
	indexes = ensure_hisat2_indexes(params.refIndexes,
									params.refGenome,
									params.refGenomeGtf,
									params.ercc,
									params.erccGtf)

	fastqs = LOAD_METASHEET(params.metaIn)

	unmapped_reads = hisat2(fastqs, indexes)

	SAVE_METASHEET(unmapped_reads, params.metaOut)
}
