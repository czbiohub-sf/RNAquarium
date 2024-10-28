include { PROCESS_UNMAPPED } from './subworkflows/process_unmapped.nf'
include { ASSEMBLE         } from './subworkflows/assembly.nf'
include { NEGATIVE_FILTER  } from './subworkflows/negative_filter.nf'
include { ALIGN_NT         } from './subworkflows/align_nt.nf'
include { ALIGN_NR         } from './subworkflows/align_nr.nf'

workflow {
    PROCESS_UNMAPPED(
        params.unmerged_accessions,
        params.bioproj_map,
        params.sra_tab_file,
        params.accession_list,
    )
    ASSEMBLE(PROCESS_UNMAPPED.out)
    NEGATIVE_FILTER(ASSEMBLE.out)
    ALIGN_NT(NEGATIVE_FILTER.out)
    ALIGN_NR(NEGATIVE_FILTER.out)
}
