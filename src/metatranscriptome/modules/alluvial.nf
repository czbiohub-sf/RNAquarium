process ALLUVIAL_PLOT {
    input:
    path nt_files, arity: "1..*"
    path nr_files, arity: "1..*"

    output:
    path 'taxonomy_hits_nonhost_alluvialplot_all.png'
    path 'taxonomy_hits_nonhost_alluvialplot_all.pdf'
    path 'taxonomy_hits_nonhost_treemap.png'
    path 'taxonomy_hits_nonhost_treemap.pdf'
    path 'taxonomy_hits_nonhost_treemap_counts.tsv'
    path 'taxonomy_hits_nonhost_fullcols_mostrecent.tsv.gz'
    path 'taxonomy_hits_nonhost_NRbetter.tsv'
    path 'taxonomy_hits_nonhost_taxonmismatch.tsv'
    path 'taxonomy_hits_viruses_fullcols_mostrecent.tsv'
    path 'taxonomy_hits_viruses0_fullcols_mostrecent.tsv'
    path 'taxonomy_hits_viruses0.tsv'
    path 'taxonomy_hits_bacteria.tsv'
    path 'taxonomy_hits_arthropoda.tsv'
    path 'taxonomy_hits_plants.tsv'
    path 'taxonomy_hits_chordates.tsv'
    path 'taxonomy_hits_fungi.tsv'
    path 'taxonomy_hits_otherEukaryota.tsv'
    path 'taxonomy_hits_SAReukaryotes.tsv'
    path 'taxonomy_hits_archaea.tsv'
    path 'taxonomy_hits_mollusca.tsv'
    path 'taxonomy_hits_annelida.tsv'
    path 'taxonomy_hits_nematoda.tsv'
    path 'taxonomy_hits_platyhelminthes.tsv'
    path 'taxonomy_hits_nonhost_list.txt'
    path 'taxonomy_hits_nonhost_alluvialplot_counts.tsv'
    path 'taxonomy_hits_nonhost_treemap_counts.tsv'
    path 'allchunks_blastn_hits_nonhost.tsv.gz'
    path 'allchunks_diamond_hits_nonhost.tsv.gz'
    path 'taxonomy_hits_nonhost_mostrecent.tsv.gz'
    path 'taxonomy_hits_viruses_list.txt'
    path 'taxonomy_hits_bacteria_list.txt'
    path 'taxonomy_hits_arthropoda_list.txt'
    path 'taxonomy_hits_plants_list.txt'
    path 'taxonomy_hits_chordates_list.txt'
    path 'taxonomy_hits_fungi_list.txt'
    path 'taxonomy_hits_otherEukaryota_list.txt'
    path 'taxonomy_hits_SAReukaryotes_list.txt'
    path 'taxonomy_hits_archaea_list.txt'
    path 'taxonomy_hits_mollusca_list.txt'
    path 'taxonomy_hits_annelida_list.txt'
    path 'taxonomy_hits_nematoda_list.txt'
    path 'taxonomy_hits_platyhelminthes_list.txt'

    script:
    """
    mkdir nt
    mkdir nr

    mv $nt_files nt/
    mv $nr_files nr/

    alluvial_treemap.r -nt nt -nr nr
    """
}
