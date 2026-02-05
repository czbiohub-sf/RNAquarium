process BLAST {
    input:
    path transcripts

    output:
    tuple path("${transcripts.simpleName}.blast.txt.gz"), path(transcripts)

    script:
    """
    OUTPUT_STAGING="\${PWD}/${transcripts.simpleName}.blast.txt.staging"
    OUTPUT="\${PWD}/${transcripts.simpleName}.blast.txt.gz"
    INPUT="\${PWD}/${transcripts}"

    # Seems like we need to be in the DB directory?
    cd /db
    blastn \
        -db /db/${params.nt_db_name} \
        -query \$INPUT \
        -taxids "${params.taxids}" \
        -outfmt "${params.outfmt}" \
        -num_threads ${task.cpus} \
        -evalue ${params.evalue} \
        -max_target_seqs "${params.max_target_seqs}" \
        -out \${OUTPUT_STAGING}
    gzip \${OUTPUT_STAGING}
    mv "\${OUTPUT_STAGING}.gz" "\${OUTPUT}"
    """
}

process CUT_BLAST_RESULTS {
    input:
    tuple path(blast_result), path(fasta)
    
    output:
    path "${blast_result.simpleName}_nonzfhum.fasta"

    script:
    """
    COLS_FILE="${blast_result.simpleName}_unique_list.txt"
    NAMES_FILE="${blast_result.simpleName}_names.txt"
    NONZFHUM_FILE="${blast_result.simpleName}_nonzfhum.txt"
    NONZFHUM_FA_FILE="${blast_result.simpleName}_nonzfhum.fasta"

    zcat $blast_result | cut -f1 | sort | uniq > \$COLS_FILE
    grep "^>" $fasta | sed 's/^>//g' > \$NAMES_FILE
    grep -v -F -f \$COLS_FILE \$NAMES_FILE > \$NONZFHUM_FILE
    seqtk subseq $fasta \$NONZFHUM_FILE >> \$NONZFHUM_FA_FILE
    """
}

process BLAST_FULL_NT {
    input:
    path transcripts

    output:
    path "${transcripts.simpleName}.blast.txt.gz"

    script:
    """
    OUTPUT_STAGING="\${PWD}/${transcripts.simpleName}.blast.txt.staging"
    OUTPUT="\${PWD}/${transcripts.simpleName}.blast.txt.gz"
    INPUT="\${PWD}/${transcripts}"

    # Seems like we need to be in the DB directory?
    cd /db
    blastn \
        -db /db/${params.ntfull_db_name} \
        -query \$INPUT \
        -outfmt "${params.full_nt_outfmt}" \
        -num_threads ${task.cpus} \
        -evalue ${params.full_nt_evalue} \
        -max_target_seqs "${params.full_nt_max_target_seqs}" \
        -out \${OUTPUT_STAGING}
    gzip \${OUTPUT_STAGING}
    mv "\${OUTPUT_STAGING}.gz" "\${OUTPUT}"
    """
}
