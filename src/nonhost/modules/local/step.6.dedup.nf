#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.publishDir = "$PWD"
params.publishIntermediate = true
params.percentLen = 100
params.minLen = 20
params.cleanupScript = ""
params.tmp = null
params.backupTmp = null
params.backupScratchHack = false
params.nxfUnstageHack = false

params.metaIn = 'step_5_sheet.csv'
params.metaOut = 'step_6_sheet.csv'
include {
	LOAD_METASHEET;
	SAVE_METASHEET;
} from './utils.nf'

process dedup {
	label 'dedup'

	input:
	tuple val(meta), path(mategz, arity: '1..2'), path(star_log)
	// when would we ever /not/ have STAR logs and want to continue this readset?

	output:
	tuple val(meta), path("Unmapped.out.mate?.filteredbyBT.dedup.fastq.gz", arity: '1..2'), emit: mates
	tuple val(meta), stdout, emit: stats

	script:
	def PERC_COEFF = params.percentLen / 100
	def PRECONDITIONS="""if ! ${task.ext.gzipCmd} -t ${mategz}; then
		>&2 echo "[Fatal error][$meta.id] At least one of the $mategz input files failed the integrity test"
	fi
	if ! [ -e $star_log ]; then
		>&2 echo "STAR final log not found: $star_log, aborting."
	fi

	read_length=\$(cat $star_log | sed '7q;d'| cut -f2)
	# check that this is a good number
	num_pat='^[0-9]+\$'
	read_length=\$(echo \$read_length | sed 's/[^0-9]*//g')
	if [ \$read_length -eq 0 ] || ! [[ \$read_length =~ \$num_pat ]] || (( \$read_length < $params.minLen ))
	then
		>&2 echo "[Fatal error][$meta.id] Cannot determine read length or read too short, got read length: \$read_length, "
		exit 1
	fi
	length_for_dedup=\$(bc <<< "scale=3; \$read_length * $PERC_COEFF")
	echo "$meta.id read length \$read_length, use \$length_for_dedup for dedup"
	"""
	def OUT_MATE1_NAME = "Unmapped.out.mate1.filteredbyBT.dedup.fastq"
	def OUT_MATE2_NAME = "Unmapped.out.mate2.filteredbyBT.dedup.fastq"
	if (!meta.single_end)
	"""
	${PRECONDITIONS}
	${task.ext.gzipCmd} -cd ${mategz[0]} > Unmapped.out.mate1.filteredbyBT.fastq
	${task.ext.gzipCmd} -cd ${mategz[1]} > Unmapped.out.mate2.filteredbyBT.fastq

	czid-dedup -i Unmapped.out.mate1.filteredbyBT.fastq -i Unmapped.out.mate2.filteredbyBT.fastq \
		-o $OUT_MATE1_NAME -o $OUT_MATE2_NAME \
		-c cluster.csv -l \$length_for_dedup

	rm cluster.csv
	${task.ext.gzipCmd} -nc $OUT_MATE1_NAME > ${OUT_MATE1_NAME}.gz.staging
	${task.ext.gzipCmd} -nc $OUT_MATE2_NAME > ${OUT_MATE2_NAME}.gz.staging
	mv ${OUT_MATE1_NAME}.gz.staging ${OUT_MATE1_NAME}.gz
	mv ${OUT_MATE2_NAME}.gz.staging ${OUT_MATE2_NAME}.gz

	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
	else if (meta.single_end)
	"""
	${PRECONDITIONS}
	${task.ext.gzipCmd} -cd ${mategz} > Unmapped.out.mate1.filteredbyBT.fastq

	czid-dedup -i Unmapped.out.mate1.filteredbyBT.fastq \
		-o $OUT_MATE1_NAME \
		-c cluster.csv -l \$length_for_dedup

	rm cluster.csv
	${task.ext.gzipCmd} -nc $OUT_MATE1_NAME > ${OUT_MATE1_NAME}.gz.staging
	mv ${OUT_MATE1_NAME}.gz.staging ${OUT_MATE1_NAME}.gz

	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
}
