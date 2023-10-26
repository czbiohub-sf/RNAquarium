#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def star_usage() {
	return """STAR index generation parameters (one of the below options is required):
Option 1:
--ref_indexes_ercc path   pre-generated STAR indexes, with ERCC spike-in controls
--ref_indexes path        pre-generated STAR indexes for reference genome, without spike-in
--genome_size n           reference genome size, in bytes. (fa file size is good enough)
Option 2:
--ref_genome path         path to ensembl reference genome fasta
	e.g. Danio_rerio.GRCz11.dna_sm.primary_assembly.fa
--ref_genome_gtf path     path to ensembl reference genome annotations
    e.g. Danio_rerio.GRCz11.108.gtf
--ercc_fa path            path to ERCC spike-in control sequences
                            (default: ERCC92.fa)
--ercc_gtf path           path to ERCC spike-in control annotations
                            (default: ERCC92.gtf)
"""
}

params.publish_dir = "$PWD"
params.publish_intermediate = true
params.genome_size = null
params.ref_indexes_ercc = null // "Danio_rerio.GRCz11.108.ERCC"
params.ref_indexes = null //"Danio_rerio.GRCz11.108"

params.ref_genome = "Danio_rerio.GRCz11.dna_sm.primary_assembly.fa"
params.ref_genome_gtf = "Danio_rerio.GRCz11.108.gtf"
params.ercc_fa = "ERCC92.fa"
params.ercc_gtf = "ERCC92.gtf"

params.star_options = ""
params.star_index_gen_options = ""

include {
	star_generate_indexes;
} from './step.0.generate_indexes.nf' params(
	publish_dir: params.publish_dir,
	star_index_gen_options: params.star_index_gen_options
)

params.meta_in = 'step_3_sheet.csv'
params.meta_out = 'step_4_sheet.csv'
include {
	LOAD_METASHEET;
	SAVE_METASHEET;
} from './utils.nf'


process star {
	label 'star'
	publishDir "$params.publish_dir/STAR_out/$meta.id", enabled: params.publish_intermediate
	
	input:
	tuple val(meta), path(fqgz)
	path indexes_dir // "Danio_rerio.GRCz11.108.ERCC"

	output:
	tuple val(meta), path("?E/Unmapped.out.mate?.gz")

	script:
	def STAR_CMD = """STAR --outFilterMultimapNmax 99999 --outFilterMismatchNmax 999 \
		--outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5 \
		--outSAMmode None \
		--clip3pNbases 0 --limitOutSJcollapsed 200000000 \
		--genomeLoad NoSharedMemory --outReadsUnmapped Fastx \
		--runThreadN ${task.cpus} \
		${params.star_options}"""
	if (!meta.single_end)
	"""
	${STAR_CMD} \
		--genomeDir ${indexes_dir} \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix PE/ \
		--readFilesIn ${fqgz[0]} ${fqgz[1]}

	gzip PE/Unmapped.out.mate*
	gzip -t PE/Unmapped.out.mate*.gz
	"""
	else if (meta.single_end)
	"""
	${STAR_CMD} \
		--genomeDir ${indexes_dir} \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix SE/ \
		--readFilesIn ${fqgz}

	gzip SE/Unmapped.out.mate1
	gzip -t SE/Unmapped.out.mate1.gz
	"""
}

def set_genome_size(n) {
	params.star_genome_size = n
}
def ensure_genome_size(explicit_size, ref_genome_fa) {
	if (explicit_size && explicit_size > 0) {
		set_genome_size(explicit_size)
		return explicit_size // user provided a genome size
	} else if (ref_genome_fa) {
		set_genome_size(ref_genome_fa.size())
		return ref_genome_fa.size() // based on fasta file size
	} else {
		log.error "--genome_size must be provided\n${star_usage()}"
		exit 1
	}
}

def ensure_star_indexes(ref_indexes, ref_indexes_ercc,
						ref_genome, ref_genome_gtf, ercc, ercc_gtf) {
	if (ref_indexes && ref_indexes_ercc) {
		(indexes, indexes2) = [file(ref_indexes_ercc), file(ref_indexes)]
	} else {
		(indexes, indexes2) = star_generate_indexes(file(ref_genome),
													file(ref_genome_gtf),
													file(ercc),
													file(ercc_gtf))
	}
	return [indexes, indexes2]
}

workflow {
	ensure_genome_size(params.genome_size, file(params.ref_genome))

	(indexes, indexes2) = ensure_star_indexes(params.ref_indexes,
											  params.ref_indexes_ercc,
											  params.ref_genome,
											  params.ref_genome_gtf,
											  params.ercc,
											  params.ercc_gtf)

	fastqs = LOAD_METASHEET(params.meta_in)

	unmapped_reads = star(fastqs, indexes)

	SAVE_METASHEET(unmapped_reads, params.meta_out)
}
