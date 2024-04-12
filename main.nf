include { merge_unmapped                       } from './modules/merge_unmapped.nf'
include { spades_single_end; spades_paired_end } from './modules/assembly.nf'
include { parse_accessions; download_sra_tab   } from './modules/accession_mapping.nf'

workflow {
    if (!params.unmerged_accessions) {
        log.error("No directory was provided for unmerged accessions!")
        exit 1
    }

    if (!params.bioproj_map) {
        if (!params.accession_list) {
            log.error("Either a bioproject mapping or accession list must be provided.")
            exit 1
        } else {
            if (!params.sra_tab_file) {
                log.info("Downloading SRA accession list...")
                sra_tab_file = download_sra_tab()
            } else {
                sra_tab_file = Channel.fromPath(params.sra_tab_file)
            }
        }

        bioproj_map = parse_accessions(
            sra_tab_file,
            Channel.fromPath(params.accession_list)
        ).first()
    } else {
        bioproj_map = Channel.fromPath(params.bioproj_map).first()
    }

    // We use .first to turn the queue channel into a value channel.
    // This is done so we can reuse the bioproj_map channel as both a collection
    // of IDs and the mapping file.
    bioproj_ids = bioproj_map.splitJson().map{ it.key }

    results = merge_unmapped(bioproj_ids, bioproj_map, params.unmerged_accessions)

    single_end_fqs = results.filter{ it[1].size() == 1 }
    paired_end_fqs = results.filter{ it[1].size() == 2 }
    both_end_fqs   = results.filter{ it[1].size() == 3 }

    // https://stackoverflow.com/a/75248731
    single_end_fqs = single_end_fqs
        .concat(
            both_end_fqs.map{ x, y -> tuple(x, y.findAll{ it =~ /PRJNA\d+_S/ }) }
        )
        .map{ id, fqs -> tuple(id, fqs[0]) }

    paired_end_fqs = paired_end_fqs
        .concat(
            both_end_fqs.map{ x, y -> tuple(x, y.findAll{ it =~ /PRJNA\d+_P/ }) }
        )
        .map{ id, fqs -> tuple(id, fqs[0], fqs[1]) }

    spades_single_end(single_end_fqs)
    spades_paired_end(paired_end_fqs)
}
