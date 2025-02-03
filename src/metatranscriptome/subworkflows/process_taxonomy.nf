include { PROCESS_BLAST } from '../modules/taxonomy.nf'

workflow PROCESS_TAXONOMY {
    take:
        blast_results

    main:
        PROCESS_BLAST(blast_results)
}
