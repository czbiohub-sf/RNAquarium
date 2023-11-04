#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.publishDir = "$PWD"
params.publishIntermediate = true
params.genomeSize = null
params.refIndexes = null

params.refGenome = "Danio_rerio.GRCz11.dna_sm.primary_assembly.fa"
params.refGenomeGtf = "Danio_rerio.GRCz11.108.gtf"
params.erccFa = "ERCC92.fa"
params.erccGtf = "ERCC92.gtf"

params.hisatOptions = ""
params.hisatIndexGenOptions = ""

include {
	hisat2_generate_indexes;
} from './step.0.generate_indexes.nf' params(
	publishDir: params.publishDir,
	hisat2IndexGenOptions: params.hisatIndexGenOptions
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
	tuple val(meta), path(fqgz)
	path indexes_dir

	output:
	tuple val(meta), path("?E/Unmapped.out.mate?.gz")

	script:
	def HISAT2_CMD = """hisat2 -q -p $task.cpus -k 1 -S /dev/null \
		-x $indexes_dir/${indexes_dir.getName()} $params.hisatOptions"""
	if (!meta.single_end)
	"""
	mkdir -p PE
	${HISAT2_CMD} \
		-1 ${fqgz[0]} -2 ${fqgz[1]} \
		--un-conc-gz PE/Unmapped.out.mate%.gz.staging
	mv PE/Unmapped.out.mate1.gz.staging PE/Unmapped.out.mate1.gz
	mv PE/Unmapped.out.mate2.gz.staging PE/Unmapped.out.mate2.gz
	"""
	else if (meta.single_end)
	"""
	mkdir -p SE
	${HISAT2_CMD} \
		-U ${fqgz} \
		--un-gz SE/Unmapped.out.mate1.gz.staging
	mv SE/Unmapped.out.mate1.gz.staging SE/Unmapped.out.mate1.gz
	"""
}

def ensure_hisat2_indexes(ref_indexes,
						  ref_genome, ref_genome_gtf, ercc, ercc_gtf) {
	if (ref_indexes) {
		indexes = file(ref_indexes)
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
