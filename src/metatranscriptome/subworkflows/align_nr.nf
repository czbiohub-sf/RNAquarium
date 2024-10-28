include { CHUNK_NONZFHUM_FASTA; DIAMOND } from '../modules/diamond.nf'

workflow ALIGN_NR {
    take:
        non_zf_hum_fa

    main:
        if (!params.nr_dir) {
            log.error("No value was provided for nr_dir!")
            exit 1
        }
        diamond_chunks = CHUNK_NONZFHUM_FASTA(
            non_zf_hum_fa.flatten(),
        )
        diamond_results = DIAMOND(diamond_chunks.flatten())

    emit:
        diamond_results
}
