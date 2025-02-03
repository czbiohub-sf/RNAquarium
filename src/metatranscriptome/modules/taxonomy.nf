process PROCESS_BLAST {
    input:
    path blast_result

    output:
    path "${blast_result.simpleName}.blastn.tab"

    script:
    """
    nt_blast_processing.r -i $blast_result
    """
}
