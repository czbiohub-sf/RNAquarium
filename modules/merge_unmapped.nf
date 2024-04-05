process merge_unmapped {
    publishDir "output/"

    input:
    val bioproject_id
    path bioproject_mapping
    path gsnap_dir

    output:
    tuple val(bioproject_id), path("${bioproject_id}_*/*.fastq.gz", arity: '1..3')

    script:
    """
    merge_unmapped.py \
        --ID $bioproject_id \
        --mapping $bioproject_mapping \
        --indir $gsnap_dir \
        --outdir .
    """
}
