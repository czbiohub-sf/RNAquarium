include { MERGE_UNMAPPED                     } from './modules/merge_unmapped.nf'
include { PARSE_ACCESSIONS; DOWNLOAD_SRA_TAB } from './modules/accession_mapping.nf'

workflow PROCESS_UNMAPPED {
    take:
        unmerged_accessions
        bioproj_map
        sra_tab_file
        accession_list
        bioproj_map

    main:
        if (!unmerged_accessions) {
            log.error("No directory was provided for unmerged accessions!")
            exit 1
        }

        if (!bioproj_map) {
            if (!accession_list) {
                log.error("Either a bioproject mapping or accession list must be provided.")
                exit 1
            }

            if (!sra_tab_file) {
                log.info("Downloading SRA accession list...")
                sra_tab_file = DOWNLOAD_SRA_TAB()
            } else {
                sra_tab_file = Channel.fromPath(sra_tab_file)
            }

            PARSE_ACCESSIONS(
                sra_tab_file,
                Channel.fromPath(accession_list)
            )
            bioproj_map = PARSE_ACCESSIONS.out.mapping.first()
        } else {
            bioproj_map = Channel.fromPath(bioproj_map).first()
        }

        // We use .first to turn the queue channel into a value channel.
        // This is done so we can reuse the bioproj_map channel as both a collection
        // of IDs and the mapping file.

        bioproj_ids = bioproj_map.splitJson().map{ it.key }
        results = MERGE_UNMAPPED(bioproj_ids, bioproj_map, unmerged_accessions)

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

    emit:
        single_end_fqs
        paired_end_fqs
}
