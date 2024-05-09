process BLAST {
    input:
    path transcripts

    output:
    path "${transcripts.simpleName}.blast.txt.gz"

    script:
    """
    OUTPUT="\${PWD}/${transcripts.simpleName}.blast.txt"
    INPUT="\${PWD}/${transcripts}"

    # Seems like we need to be in the DB directory?
    cd /db
    blastn \
        -db /db/nt \
        -query \$INPUT \
        -taxids "${params.taxids}" \
        -outfmt "${params.outfmt}" \
        -num_threads ${task.cpus} \
        -evalue ${params.evalue} \
        -max_target_seqs "${params.max_target_seqs}" \
        -out \$OUTPUT
    gzip \$OUTPUT
    """
}

process CONCAT_BLAST {
    input:
    path blast_results, arity: '1..*'

    output:
    path "concat_blast.txt.gz"

    script:
    """
    # TODO: Use params.outfmt to create header.
    #       Need to account for different outfmts (i.e. tab vs CSV)
    #       Outfmt 7 also adds comment lines to each file so maybe don't concat?
    cat ${blast_results} > concat_blast.txt.gz
    """
}
