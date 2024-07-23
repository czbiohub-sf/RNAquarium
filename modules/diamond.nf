process CHUNK_NONZFHUM_FASTA {
    input:
    path transcript

    output:
    path "${transcript.simpleName}.*.fasta", arity: params.diamond_chunks

    script:
    num_chunks = params.diamond_chunks
    """
    echo ${transcript.simpleName}
    alphabet="abcdefghijklmnopqrstuvwxyz"
    seqkit split2 -p ${num_chunks} ${transcript }
    file_chunks=(\$(ls ${transcript}.split))

    i=1
    for f in ${transcript}.split/*.fasta
    do
        echo \$f
        letter=\${alphabet:i-1:1}
        mv \$f ${transcript.simpleName}.\$letter.fasta
        i=\$((i+1))
    done

    echo \$(ls)
    """
}

process DIAMOND {
    input:
    path transcripts

    output:
    path "${transcripts.simpleName}.diamond.txt.gz"

    script:
    """
    OUTPUT="\${PWD}/${transcripts.simpleName}.diamond.txt"
    INPUT="\${PWD}/${transcripts}"
    cd /db

    diamond blastx \
        --ultra-sensitive \
        --db /db/nr \
        --threads $task.cpus \
        --query \$INPUT \
        --outfmt 6 qseqid sseqid staxids sscinames sskingdoms pident length mismatch qcovhsp gapopen qstart qend sstart send evalue bitscore \
        --top 3 \
        --evalue 0.05 \
        --out \$OUTPUT
    gzip \$OUTPUT
    """
}

process CHUNK_FULL_NT {
    input:
    path transcripts, arity: '1..*'

    output:
    path "full_nt.chunk_*.fasta"

    script:
    num_chunks = params.full_nt_blast_chunks
    """
    """
}
