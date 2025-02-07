include { PROCESS_BLAST; PROCESS_DIAMOND } from '../modules/taxonomy.nf'
include { ALLUVIAL_PLOT                  } from '../modules/alluvial.nf'

workflow PROCESS_TAXONOMY {
    take:
        blast_results
        diamond_results

    main:
        blast_tax = PROCESS_BLAST(blast_results.filter{ it.size() > 100 }) | collect
        diamond_tax = PROCESS_DIAMOND(diamond_results.filter{ it.size() > 100 }) | collect
        ALLUVIAL_PLOT(blast_tax, diamond_tax)
}
