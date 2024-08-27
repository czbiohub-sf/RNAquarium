process CHUNK_NONZFHUM_FASTA {
    input:
    path transcript
    val diamond_chunks

    output:
    path "${transcript.simpleName}.*.fasta", arity: params.diamond_chunks

    script:
    """
    echo ${transcript.simpleName}
    alphabet="abcdefghijklmnopqrstuvwxyz"
    seqkit split2 -p ${diamond_chunks} ${transcript }
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

// Use baseName over simpleName b/c transcripts are named chunk_XXX.Y.fasta
//     where Y is a letter corresponding to the sub-chunk (e.g. a-j for 1-10)
process DIAMOND {
    input:
    path transcripts

    scratch true

    output:
    path "${transcripts.baseName}.diamond.txt.gz"

    script:
    """
    OUTPUT_STAGING="\${PWD}/${transcripts.baseName}.diamond.txt.staging"
    OUTPUT="\${PWD}/${transcripts.baseName}.diamond.txt.gz"
    INPUT="\${PWD}/${transcripts}"

    diamond blastx \
        --ultra-sensitive \
        --db /db/nr \
        --threads $task.cpus \
        --query \$INPUT \
        --outfmt 6 qseqid sseqid staxids sscinames sskingdoms pident length mismatch qcovhsp gapopen qstart qend sstart send evalue bitscore stitle \
        --top 3 \
        --evalue 0.05 \
        --out \${OUTPUT_STAGING}
    gzip \${OUTPUT_STAGING}
    mv "\${OUTPUT_STAGING}.gz" "\${OUTPUT}"
    """
}
