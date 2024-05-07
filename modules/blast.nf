process BLAST {
    tag "${transcripts.simpleName}"
    container 'staphb/blast:2.15.0'
    containerOptions "--mount type=bind,src=${params.nt_dir},dst=/db"

    input:
    path transcripts

    output:
    path "${transcripts.simpleName}.blast.txt"

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
    """
}
