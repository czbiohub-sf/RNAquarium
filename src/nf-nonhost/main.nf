#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.accessions_list = "SRA_accession_list.test.txt"

include {
	prefetch;
	fastq_dump;
} from './step.1.nf'
include {
	filter_barcodes;
	fastp;
	priceseqfilter;
} from './step.2.nf'

workflow {
	main:
	// step 1
	accessions = channel.fromPath(params.accessions_list).splitText()
	fastqs = prefetch(accessions) | fastq_dump

	// step 2
	filter_barcodes(fastqs) | fastp | priceseqfilter
	
}