include { CHUNK_ASSEMBLED_FASTAS   } from '../modules/preblast.nf'
include { BLAST; CUT_BLAST_RESULTS } from '../modules/blast.nf'

workflow NEGATIVE_FILTER {
    take:
        assembled_transcripts

    main:
        if (!params.nt_dir) {
            log.error("No value was provided for nt_dir!")
            exit 1
        }
        if (!params.nt_db_name) {
            log.error("No value was provided for nt_db_name!")
            exit 1
        }
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
