include { MERGE_UNMAPPED                          } from './modules/merge_unmapped.nf'
include { SPADES_SINGLE_END; SPADES_PAIRED_END    } from './modules/assembly.nf'
include { PARSE_ACCESSIONS; DOWNLOAD_SRA_TAB      } from './modules/accession_mapping.nf'
include { CHUNK_ASSEMBLED_FASTAS2                 } from './modules/preblast.nf'
include { BLAST; CUT_BLAST_RESULTS; BLAST_FULL_NT } from './modules/blast.nf'
include { CHUNK_NONZFHUM_FASTA } from './modules/diamond.nf'

workflow {
    if (!params.unmerged_accessions) {
        log.error("No directory was provided for unmerged accessions!")
        exit 1
    }

    if (!params.bioproj_map) {
        if (!params.accession_list) {
            log.error("Either a bioproject mapping or accession list must be provided.")
            exit 1
        }

        if (!params.sra_tab_file) {
            log.info("Downloading SRA accession list...")
            sra_tab_file = DOWNLOAD_SRA_TAB()
        } else {
            sra_tab_file = Channel.fromPath(params.sra_tab_file)
        }

        PARSE_ACCESSIONS(
            sra_tab_file,
            Channel.fromPath(params.accession_list)
        )
        bioproj_map = PARSE_ACCESSIONS.out.mapping.first()
    } else {
        bioproj_map = Channel.fromPath(params.bioproj_map).first()
    }

    // We use .first to turn the queue channel into a value channel.
    // This is done so we can reuse the bioproj_map channel as both a collection
    // of IDs and the mapping file.
    bioproj_ids = bioproj_map.splitJson().map{ it.key }

    results = MERGE_UNMAPPED(bioproj_ids, bioproj_map, params.unmerged_accessions)

    single_end_fqs = results.filter{ it[1].size() == 1 }
    paired_end_fqs = results.filter{ it[1].size() == 2 }
    both_end_fqs   = results.filter{ it[1].size() == 3 }

    // https://stackoverflow.com/a/75248731
    single_end_fqs = single_end_fqs
        .concat(
            both_end_fqs.map{ x, y -> tuple(x, y.findAll{ it =~ /PRJ[A-Z]{2}\d+_S/ }) }
        )
        .map{ id, fqs -> tuple(id, fqs[0]) }

    paired_end_fqs = paired_end_fqs
        .concat(
            both_end_fqs.map{ x, y -> tuple(x, y.findAll{ it =~ /PRJ[A-Z]{2}\d+_P/ }) }
        )
        .map{ id, fqs -> tuple(id, fqs[0], fqs[1]) }

    single_end_transcripts = SPADES_SINGLE_END(single_end_fqs)
    paired_end_transcripts = SPADES_PAIRED_END(paired_end_fqs)

    all_transcripts = single_end_transcripts
        .mix(paired_end_transcripts)
        .collect()

    // chunked_transcripts = CHUNK_ASSEMBLED_FASTAS(all_transcripts).flatten()
    chunked_transcripts = CHUNK_ASSEMBLED_FASTAS2(all_transcripts).flatten()
    blast_results = BLAST(chunked_transcripts)
    non_zf_hum_fa = CUT_BLAST_RESULTS(blast_results)
    full_nt_blast_results = BLAST_FULL_NT(non_zf_hum_fa)

    diamond_chunks = non_zf_hum_fa.flatten() | CHUNK_NONZFHUM_FASTA
}
