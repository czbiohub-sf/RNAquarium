process BLAST {
    tag "${transcripts}"
    container 'staphb/blast:2.15.0'
    // TODO: Update to make mount arg agnostic to containerization software
    containerOptions "--bind ${params.nt_dir}:/db"

    input:
    path transcripts

    output:
    path "${transcripts}.blast.txt"

    script:
    """
    OUTPUT="\${PWD}/${transcripts}.blast.txt"
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
