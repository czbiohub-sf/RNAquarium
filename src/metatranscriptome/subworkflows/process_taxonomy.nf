include { PROCESS_BLAST; PROCESS_DIAMOND } from '../modules/taxonomy.nf'

workflow PROCESS_TAXONOMY {
    take:
        blast_results
        diamond_results

    main:
        // PROCESS_BLAST(blast_results)
        PROCESS_DIAMOND(diamond_results)
}
