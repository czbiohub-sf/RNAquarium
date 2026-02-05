process SPADES_SINGLE_END {
    tag "${bioproj_id}"
    label 'use_scratch'

    input:
    tuple val(bioproj_id), path(fq)

    output:
    path "${bioproj_id}_S.transcripts.fasta"

    script:
    """
    spades.py --rna -s $fq -o .
    if [ -f transcripts.fasta ]; then
        sed -i 's/^>NODE/>${bioproj_id}_S_NODE/g' transcripts.fasta
        mv transcripts.fasta ${bioproj_id}_S.transcripts.fasta
    else
        echo 'No transcripts.fasta file found!'
        exit 100
    fi
    """
}

process SPADES_PAIRED_END {
    tag "${bioproj_id}"
    label 'use_scratch'

    input:
    tuple val(bioproj_id), path(fq1), path(fq2)

    output:
    path "${bioproj_id}_P.transcripts.fasta"

    script:
    """
    spades.py --rna -1 $fq1 -2 $fq2 -o .
    if [ -f transcripts.fasta ]; then
        sed -i 's/^>NODE/>${bioproj_id}_P_NODE/g' transcripts.fasta
        mv transcripts.fasta ${bioproj_id}_P.transcripts.fasta
    else
        echo 'No transcripts.fasta file found!'
        exit 100
    fi
    """
}
