include { CHUNK_NONZFHUM_FASTA; DIAMOND           } from './modules/diamond.nf'

workflow {
    take:
        non_zf_hum_fa

    main:
        diamond_chunks = non_zf_hum_fa.flatten() | CHUNK_NONZFHUM_FASTA
        diamond_results = diamond_chunks | DIAMOND

    emit:
        diamond_results
}
