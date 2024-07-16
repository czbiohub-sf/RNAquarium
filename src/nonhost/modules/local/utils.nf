nextflow.enable.dsl=2

params.cleanupScript = ""
params.tmp = null
params.backupTmp = null
params.backupScratchHack = false


// we may want to validate that mate fq/bam/etc. filepaths exist here
// (as the error may be more legible than whatever a process throws
workflow LOAD_METASHEET {
	take:
	metasheet_name

	main:
		Channel.fromPath(metasheet_name)
		.splitCsv(header: true)
		.map { row ->
			meta = [id: row.id,
					reads: row.reads.toInteger(),
					sra_size: row.sra_size.toInteger(),
					readlen: row.readlen.toInteger(),
					fastq_size: row.fastq_size.toInteger(),
					single_end: row.single_end ]
			if (row.mate_2 > "")
				mate_paths = [row.mate_1, row.mate_2]
			else
				mate_paths = [row.mate_1]
			[ meta, mate_paths ]
		}
		.set{mates}

	emit:
	mates
}

workflow SAVE_METASHEET {
	take:
	mates_tuple
	metasheet_name

	main:
	metasheet = file(metasheet_name)
	metasheet.text = "id,reads,sra_size,readlen,fastq_size,single_end,mate_1,mate_2"
	mates_tuple.each { meta, mates ->
		metasheet << "${meta.id},${meta.reads},${meta.sra_size},"
		metasheet << "${meta.readlen},${meta.fastq_size},${meta.single_end},"
		if (meta.single_end)
			metasheet << "${mates},\n"
		else
			metasheet << "${mates[0]},${mates[1]}\n"
	}

	emit:
	metasheet_name
}


process cleanup_branched {
	input:
	tuple val(meta), path(branchin1), path(branchin2)

	script:
	"""
	cleanup="${meta.cleanup}"
	${params.cleanupScript}
	"""
}