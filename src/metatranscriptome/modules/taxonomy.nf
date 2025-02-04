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

process PROCESS_DIAMOND {
    executor 'local'

    input:
    // Example input: chunk_010_nonzfhum.i.diamond.txt.gz
    path diamond_result

    output:
    path { renameFile(diamond_result) }

    script:
    """
    nr_diamond_processing.r -i $diamond_result
    """
}

// Helper function to rename diamond result while retaining subchunk letter
def renameFile(file) {
    def pattern = /^(.+_nonzfhum)\.([a-z])\.diamond\.txt\.gz$/
    def match = file.getName() =~ pattern
    if (!match) {
        error "Invalid filename format: ${file.getName()}"
    }
    return "${match[0][1]}.${match[0][2]}.diamond_hits.tab"
}
