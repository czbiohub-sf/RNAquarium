// Look into parallelizing instead of using a for loop in bash
// Maybe create a process that takes in num_chunks, transcripts, chunk_num
//     and calculate the correct transcript splitting for each chunk separately
process CHUNK_ASSEMBLED_FASTAS {
    input:
    path transcripts, arity: '1..*'

    output:
    path "chunk_*.fasta"

    script:
    num_transcripts = transcripts.size()
    num_chunks = params.blast_chunks
    """
    TRANSCRIPTS=(${transcripts})
    BASE_CHUNK_SIZE=\$(( ${num_transcripts} / ${num_chunks} ))
    EXTRA_TRANSCRIPTS=\$(( ${num_transcripts} % ${num_chunks} ))

    START=0

    for (( i=1; i<=${num_chunks}; i++ )); do
        # Add 'remainder' transcripts sequentially starting from the first chunk
        if [[ \$i -le \$EXTRA_TRANSCRIPTS ]]; then
            CHUNK_SIZE=\$(( \$BASE_CHUNK_SIZE + 1 ))
        else
            CHUNK_SIZE=\$BASE_CHUNK_SIZE
        fi

        SUB_TRANSCRIPTS=("\${TRANSCRIPTS[@]:\$START:\$CHUNK_SIZE}")

        if [[ \${#SUB_TRANSCRIPTS[@]} -ne 0 ]]; then
            cat \$SUB_TRANSCRIPTS | seqkit shuffle - > "chunk_\$i.fasta"
        fi

        START=\$(( START + CHUNK_SIZE ))
    done
    """
}

process CHUNK_ASSEMBLED_FASTAS2 {
    input:
    path transcripts, arity: '1..*'

    output:
    path "chunk_*.fasta"

    script:
    num_chunks = params.blast_chunks
    """
    for transcript in ${transcripts}
    do
        seqkit split2 -p ${num_chunks} \$transcript 
    done

    for chunk in {001..${num_chunks}}
    do
        cat *.split/*part_\$chunk.fasta > chunk_\$chunk.fasta
    done
    """
}
