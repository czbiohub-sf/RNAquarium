process merge_unmapped {
    publishDir "${params.publish_dir}/merged_nonhost/"

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
