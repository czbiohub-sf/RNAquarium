#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.publishDir = "$PWD"
params.publishIntermediate = true
params.genomeSize = null
params.maxMismatch = 0.3
params.cleanupScript = ""
params.tmp = null
params.backupTmp = null
params.backupScratchHack = false
params.nxfUnstageHack = false

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
)

process snap {
	label 'snap'

	input:
	tuple val(meta), path(mategz, arity: '1..2')
	path index_dir

	output:
	tuple val(meta), path("snap_out.sam"), emit: sam

	script:
	def index_name = file(index_dir).getName()
	if (!meta.single_end)
	"""
	snap-aligner paired ${index_dir} ${mategz[0]} ${mategz[1]} \
		-t ${task.cpus} -o snap_out.staging.sam  -x -f -xf 2.0x
	mv snap_out.staging.sam snap_out.sam
	"""
	else if (meta.single_end)
	"""
	snap-aligner single ${index_dir} ${mategz} \
		-t ${task.cpus} -o snap_out.staging.sam -x -f -xf 2.0
	mv snap_out.staging.sam snap_out.sam
	"""
}

process gsnap {
	label 'gmap'

	def ALIGNER = "gsnap"
	def SUFFIX_OK = "filteredbyBT.dedup.gsnapFiltered.fastq"
	def SUFFIX_NG = "filteredbyBT.dedup.gsnapSkipped.fastq"
	def SAM_NAME = "${ALIGNER}_out.sam"

	input:
	tuple val(meta), path("m?.fq.gz", arity: '1..2')
	path index_dir

	output:
	tuple val(meta), path("Unmapped.out.mate?.*.gz", arity: '1..2'), emit: mates
	tuple val(meta), path("${SAM_NAME}"), emit: sam
	tuple val(meta), path("${ALIGNER}.stats.txt"), emit: stats

	script:
	def gsnap_gmap_bin = params.genomeSize < 2**32 ? "gsnap" : "gmap"
	def index_name = file(index_dir).getName()
	def ALIGNER_CMD = """${gsnap_gmap_bin} -A sam \
		--batch=4 \
		--use-shared-memory=0 \
		--maxsearch=128 \
		--npaths=1 `# maximum paths to print` \
		--ordered -t ${task.cpus} \
		--max-mismatches ${params.maxMismatch} \
		-D . -d ${index_name} \
		-o ${ALIGNER}.staging.sam """
	def ALIGNER_CMD_PE = """${ALIGNER_CMD} m1.fq m2.fq"""
	def ALIGNER_CMD_SE = """${ALIGNER_CMD} m1.fq"""

	// filter settings
	def cond = params.retainMixed ?
		(meta.single_end ? '!flag.secondary && flag.unmap' : '!flag.secondary && (flag.unmap || flag.munmap)') :
		(meta.single_end ? '!flag.secondary && flag.unmap' : '!flag.secondary && (flag.unmap && flag.munmap)')
	def NAMES = "${ALIGNER}_unmapped_names.txt"
	def SAMSTATS_CMD = """samtools view -@ ${task.cpus} ${SAM_NAME} | cut -f2 | sort | uniq -c > ${ALIGNER}.stats.txt"""
	def GET_NAMES_CMD = """samtools view -@ ${task.cpus} -e '${cond}' ${SAM_NAME} | cut -f1  > ${NAMES}"""
	def FILTER_CMD = """LC_ALL=C fastq-namefilter ${NAMES} -"""

	"""
	${!meta.single_end
	? "${task.ext.gzipCmd} -kcd m1.fq.gz > m1.fq ; ${task.ext.gzipCmd} -kcd m2.fq.gz > m2.fq"
	: "${task.ext.gzipCmd} -kcd m1.fq.gz > m1.fq"
	}
	set +e  # suppress terminate-on-error
	${!meta.single_end ? ALIGNER_CMD_PE : ALIGNER_CMD_SE}
	set -e  # resume terminate on error, check error and clear outfile.
	if [[ \$? > 0 ]]; then  # gsnap failed, pass through reads unchanged...
		:> ${SAM_NAME}
		:> ${ALIGNER}.stats.txt
		for i in ${!meta.single_end ? "{1..2}" : "1"} ; do
			cp m\${i}.fq.gz Unmapped.out.mate\${i}.${SUFFIX_NG}.gz.staging
		done
	else
		mv ${ALIGNER}.staging.sam ${SAM_NAME}
		${SAMSTATS_CMD}
		${GET_NAMES_CMD}
		for i in ${!meta.single_end ? "{1..2}" : "1"} ; do
			<m\${i}.fq ${FILTER_CMD} | ${task.ext.gzipCmd} -nc > Unmapped.out.mate\${i}.${SUFFIX_OK}.gz.staging
		done
	fi
	for file in *.gz.staging ; do mv \$file \${file/.gz.staging/.gz} ; done

	cleanup="${meta.cleanup}"
	${params.cleanupScript}
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

def ensure_snap_indexes(ref_indexes,
						ref_genome, ercc) {
	if (ref_indexes
		&& (indexes = file(ref_indexes))
		&& indexes.exists()) {
	} else {
		indexes = snap_generate_indexes(file(ref_genome),
										file(ercc))
	}
	return indexes
}
