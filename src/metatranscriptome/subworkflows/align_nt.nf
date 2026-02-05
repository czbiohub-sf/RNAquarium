include { BLAST_FULL_NT } from '../modules/blast.nf'

workflow ALIGN_NT {
    take:
        non_zf_hum_fa

    main:
        if (!params.ntfull_dir) {
            log.error("No value was provided for ntfull_dir!")
            exit 1
        }
        if (!params.ntfull_db_name) {
            log.error("No value was provided for ntfull_db_name!")
            exit 1
        }
        full_nt_blast_results = BLAST_FULL_NT(non_zf_hum_fa)

    emit:
        full_nt_blast_results
}
