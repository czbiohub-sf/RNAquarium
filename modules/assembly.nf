process spades_single_end {
    publishDir 'single_end', mode: 'copy'
    container 'docker://staphb/spades'

    memory '32 GB'

    input:
    tuple val(bioproj_id), path(fq)

    output:
    path "${bioproj_id}_S.transcripts.fasta"

    script:
    """
    spades.py --rna -s $fq -o .
    mv transcripts.fasta ${bioproj_id}_S.transcripts.fasta
    """
}

process spades_paired_end {
    publishDir 'paired_end', mode: 'copy'
    container 'docker://staphb/spades'

    memory '32 GB'

    input:
    tuple val(bioproj_id), path(fq1), path(fq2)

    output:
    path "${bioproj_id}_P.transcripts.fasta"

    script:
    """
    spades.py --rna -1 $fq1 -2 $fq2 -o .
    mv transcripts.fasta ${bioproj_id}_P.transcripts.fasta
    """
}
