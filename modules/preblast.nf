process CHUNK_ASSEMBLED_FASTAS {
    input:
    path transcripts, arity: '1..*'
    val blast_chunks

    output:
    path "chunk_*.fasta"

    script:
    """
    for transcript in ${transcripts}
    do
        seqkit split2 -p ${blast_chunks} \$transcript 
    done

    for chunk in {001..${blast_chunks}}
    do
        cat *.split/*part_\$chunk.fasta > chunk_\$chunk.fasta
    done
    """
}
