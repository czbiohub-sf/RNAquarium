nextflow.enable.dsl=2
import nextflow.util.SysHelper


params.publishDir = "$PWD"
params.publishIntermediate = true
// https://ftp.ensembl.org/pub/release-108/fasta/danio_rerio/dna/
// https://ftp.ensembl.org/pub/release-108/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz
params.refGenome = "Danio_rerio.GRCz11.dna_sm.primary_assembly.fa"
// https://ftp.ensembl.org/pub/release-108/gtf/danio_rerio/
// https://ftp.ensembl.org/pub/release-108/gtf/danio_rerio/Danio_rerio.GRCz11.108.gtf.gz
params.refGenomeGtf = "Danio_rerio.GRCz11.108.gtf"
// https://tools.thermofisher.com/content/sfs/manuals/ERCC92.zip
params.erccFa = "ERCC92.fa"
params.erccGtf = "ERCC92.gtf"

params.starSjdbOverhang = 100

params.hisatUseTranscript = false

params.starIndexGenOptions = ""
params.hisatIndexGenOptions = ""

process star_generate_indexes {
	label 'star'
	cache true
	
	input:
	path ref_genome_fa
	path ref_genome_gtf
	path ERCC_fa
	path ERCC_gtf

	output:
	path("*_indexes.ERCC"), emit: ercc_indexes
	path("*_indexes"), emit: indexes

	script:
	def genome_name = ref_genome_fa.getSimpleName()
	def STAR_INDEXGEN_CMD = """STAR --runMode genomeGenerate --runThreadN ${task.cpus} \
		--sjdbOverhang ${params.starSjdbOverhang} $params.starIndexGenOptions """
	"""
	cat $ref_genome_fa $ERCC_fa > dna_sm.primary_assembly_ERCC.fa
	cat $ref_genome_gtf $ERCC_gtf > indexes_ERCC.gtf

	${STAR_INDEXGEN_CMD} --genomeFastaFiles $ref_genome_fa --sjdbGTFfile $ref_genome_gtf \
		--genomeDir	./${genome_name}_indexes

	${STAR_INDEXGEN_CMD} --genomeFastaFiles dna_sm.primary_assembly_ERCC.fa \
		--sjdbGTFfile indexes_ERCC.gtf --genomeDir ./${genome_name}_indexes.ERCC
	"""
}

process hisat2_generate_indexes {
	label 'hisat2'
	cache true

	input:
	path ref_genome_fa
	path ref_genome_gtf
	path ERCC_fa
	path ERCC_gtf

	output:
	path("*_genome")

	// untested
	script:
	def genome_name = ref_genome_fa.getSimpleName()
	if (params.hisatUseTranscript && ref_genome_gtf.exists())
	"""
	# can potentially be done in parallel
	cat $ref_genome_fa $ERCC_fa > dna_sm.primary_assembly_ERCC.fa
	cat $ref_genome_gtf $ERCC_gtf > indexes_ERCC.gtf

	#hisat2_extract_splice_sites.py indexes_ERCC.gtf > genome.ss
	hisat2_extract_exons.py indexes_ERCC.gtf > genome.exon

	mkdir -p ${genome_name}_genome
	hisat2-build -q -p $task.cpus $params.hisatIndexGenOptions \
		--exon genome.exon \
		dna_sm.primary_assembly_ERCC.fa ${genome_name}_genome/${genome_name}_genome
	"""
	else
	"""
	echo '[hisat2_generate_indexes] building index without GTF exon/splice graph'
	cat $ref_genome_fa $ERCC_fa > dna_sm.primary_assembly_ERCC.fa

	mkdir -p ${genome_name}_genome
	hisat2-build -q -p $task.cpus $params.hisatIndexGenOptions \
		dna_sm.primary_assembly_ERCC.fa ${genome_name}_genome/${genome_name}_genome
	"""
}

workflow {
	(star_index, star_index2) = star_generate_indexes(file(params.refGenome),
													  file(params.refGenomeGtf),
													  file(params.ercc),
													  file(params.erccGtf))
	
	hisat_index = hisat2_generate_indexes(file(params.refGenome),
										  file(params.refGenomeGtf),
										  file(params.ercc),
										  file(params.erccGtf))
}