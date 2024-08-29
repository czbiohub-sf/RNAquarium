include { CHUNK_ASSEMBLED_FASTAS   } from '../modules/preblast.nf'
include { BLAST; CUT_BLAST_RESULTS } from '../modules/blast.nf'

workflow NEGATIVE_FILTER {
    take:
        assembled_transcripts

    main:
        chunked_transcripts = (
            CHUNK_ASSEMBLED_FASTAS(
                assembled_transcripts,
            ).flatten()
        )
        blast_results = BLAST(chunked_transcripts)
        non_zf_hum_fa = CUT_BLAST_RESULTS(blast_results)

    emit:
        non_zf_hum_fa
}
