include { BLAST_FULL_NT } from '../modules/blast.nf'

workflow ALIGN_NT {
    take:
        non_zf_hum_fa

    main:
        full_nt_blast_results = BLAST_FULL_NT(non_zf_hum_fa)

    emit:
        full_nt_blast_results
}
