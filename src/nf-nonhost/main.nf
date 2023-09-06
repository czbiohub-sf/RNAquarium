#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.accessions_list = "SRA_accession_list.test.txt"

include {
	prefetch;
	fastq_dump;
} from './step.1.nf'
include {
	filter_barcodes;
	fastp_pe;
	fastp_se;
	//price
} from './step.2.nf'

workflow {
	main:
	// step 1
	accessions = channel.fromPath(params.accessions_list).splitText()
	fastqs = prefetch(accessions) | fastq_dump

	// step 2
	filter_barcodes(fastqs)
		.branch {
			SE: it.size() != 2
			PE: it.size() == 2
		}
		.set { fastqs_filtered }
	fastp_pe(fastqs_filtered.PE)
	fastp_se(fastqs_filtered.SE)

	
}
	