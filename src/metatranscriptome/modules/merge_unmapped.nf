process MERGE_UNMAPPED {
    tag "${bioproject_id}"
    publishDir "${params.publish_dir}/merged_nonhost/"

    container 'docker://python:3.10'
    conda 'python=3.10'

    debug true
    errorStrategy 'ignore'
    cache 'lenient'

    input:
    val bioproject_id
    path bioproject_mapping
    path unmerged_accessions

    output:
    tuple val(bioproject_id), path("${bioproject_id}_*/*.fastq.gz", arity: '1..3')

    script:
    """
    merge_unmapped.py \
        --ID $bioproject_id \
        --mapping $bioproject_mapping \
        --indir $unmerged_accessions \
        --outdir .
    """
}
