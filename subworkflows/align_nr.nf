include { CHUNK_NONZFHUM_FASTA; DIAMOND } from '../modules/diamond.nf'

workflow ALIGN_NR {
    take:
        non_zf_hum_fa

    main:
        diamond_chunks = CHUNK_NONZFHUM_FASTA(
            non_zf_hum_fa.flatten(),
        )
        diamond_results = DIAMOND(diamond_chunks.flatten())

    emit:
        diamond_results
}
