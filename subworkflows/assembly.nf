include { SPADES_SINGLE_END; SPADES_PAIRED_END } from './modules/assembly.nf'

workflow SPADES {
    take:
        single_end_fqs
        paired_end_fqs

    main:
        single_end_transcripts = SPADES_SINGLE_END(single_end_fqs)
        paired_end_transcripts = SPADES_PAIRED_END(paired_end_fqs)
        all_transcripts = single_end_transcripts
            .mix(paired_end_transcripts)
            .collect()

    emit:
        assembled_transcripts = all_transcripts
}
