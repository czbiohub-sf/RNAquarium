process SPADES_SINGLE_END {
    tag "${bioproj_id}"
    publishDir "${params.publish_dir}/single_end"

    container 'docker://staphb/spades'
    conda 'bioconda::spades=3.15.5'

    memory { 16.GB * task.attempt }
    maxRetries 3
    errorStrategy { task.exitStatus == 9 ? 'ignore' : 'retry' }
    scratch true

    input:
    tuple val(bioproj_id), path(fq)

    output:
    path "${bioproj_id}_S.transcripts.fasta"

    script:
    """
    spades.py --rna -s $fq -o .
    if [ -f transcripts.fasta ]; then
        sed -i 's/>NODE/>${bioproj_id}_S_NODE/g' transcripts.fasta
        mv transcripts.fasta ${bioproj_id}_S.transcripts.fasta
    else
        echo 'No transcripts.fasta file found!'
        exit 9
    fi
    """
}

process SPADES_PAIRED_END {
    tag "${bioproj_id}"
    publishDir "${params.publish_dir}/paired_end"

    container 'docker://staphb/spades'
    conda 'bioconda::spades=3.15.5'

    memory { 32.GB * task.attempt }
    maxRetries 3
    errorStrategy { task.exitStatus == 9 ? 'ignore' : 'retry' }
    scratch true

    input:
    tuple val(bioproj_id), path(fq1), path(fq2)

    output:
    path "${bioproj_id}_P.transcripts.fasta"

    script:
    """
    spades.py --rna -1 $fq1 -2 $fq2 -o .
    if [ -f transcripts.fasta ]; then
        sed -i 's/>NODE/>${bioproj_id}_P_NODE/g' transcripts.fasta
        mv transcripts.fasta ${bioproj_id}_P.transcripts.fasta
    else
        echo 'No transcripts.fasta file found!'
        exit 9
    fi
    """
}
