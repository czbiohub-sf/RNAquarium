#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.publish_dir = "$PWD"
params.publish_intermediate = true
params.genome_size = null
params.ref_indexes = null

params.ref_genome = "Danio_rerio.GRCz11.dna_sm.primary_assembly.fa"
params.ref_genome_gtf = "Danio_rerio.GRCz11.108.gtf"
params.ercc_fa = "ERCC92.fa"
params.ercc_gtf = "ERCC92.gtf"

params.hisat2_options = ""
params.hisat2_index_gen_options = ""

include {
	hisat2_generate_indexes;
} from './step.0.generate_indexes.nf' params(
	publish_dir: params.publish_dir,
	hisat2_index_gen_options: params.hisat2_index_gen_options
)

params.meta_in = 'step_2_sheet.csv'
params.meta_out = 'step_3_sheet.csv'
include {
	LOAD_METASHEET;
	SAVE_METASHEET;
} from './utils.nf'

process hisat2 {
	label 'hisat2'
	publishDir "$params.publish_dir/hisat2_out/$meta.id", enabled: params.publish_intermediate

	input:
	tuple val(meta), path(fqgz)
	path indexes_dir

	output:
	tuple val(meta), path("?E/Unmapped.out.mate?.gz")

	script:
	def HISAT2_CMD = """hisat2 -q -p $task.cpus -k 1 -S /dev/null \
		-x $indexes_dir/${indexes_dir.getName()} $params.hisat2_options"""
	if (!meta.single_end)
	"""
	mkdir -p PE
	${HISAT2_CMD} \
		-1 ${fqgz[0]} -2 ${fqgz[1]} \
		--un-conc-gz PE/Unmapped.out.mate%.gz
	"""
	else if (meta.single_end)
	"""
	mkdir -p SE
	${HISAT2_CMD} \
		-U ${fqgz} \
		--un-gz SE/Unmapped.out.mate1.gz
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
	indexes = ensure_hisat2_indexes(params.ref_indexes,
									params.ref_genome,
									params.ref_genome_gtf,
									params.ercc,
									params.ercc_gtf)

	fastqs = LOAD_METASHEET(params.meta_in)

	unmapped_reads = hisat2(fastqs, indexes)

	SAVE_METASHEET(unmapped_reads, params.meta_out)
}
