process SPADES_SINGLE_END {
    publishDir "${params.publish_dir}/single_end"
    container 'docker://staphb/spades'
    conda 'bioconda::spades=3.15.5'

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

process SPADES_PAIRED_END {
    publishDir "${params.publish_dir}/paired_end"
    container 'docker://staphb/spades'
    conda 'bioconda::spades=3.15.5'

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
