process PARSE_ACCESSIONS {
    input:
    path tabfile
    path accfile

    output:
    path "bioproject_accession_mapping.json", emit: mapping
    path "unmapped_accessions.txt", emit: unmapped_accs

    script:
    """
    parse_accessions.py \
        --tab $tabfile \
        --acc $accfile \
        --outfile bioproject_accession_mapping.json \
        --unmapped unmapped_accessions.txt
    """
}

process DOWNLOAD_SRA_TAB {
    output:
    path "SRA_Accessions.tab"

    script:
    """
    wget ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab
    """
}
