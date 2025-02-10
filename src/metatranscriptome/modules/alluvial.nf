process ALLUVIAL_PLOT {
    input:
    path nt_files, arity: "1..*"
    path nr_files, arity: "1..*"

    output:
    path "allchunks_blastnclustered_hits_nonzfhum.tsv"
    path "allchunks_diamond_hits_nonzfhum.tsv"
    path "allchunks_blastnanddiamond_hits_nonzfhum_*.tsv"
    path "allchunks_blastnanddiamond_hits_*_list.txt"
    path "allchunks_blastnanddiamond_hits_*.tsv"
    path "allchunks_alluvialplot_counts.tsv"
    path "RNaquarium_allchunks_blastnanddiamond_allhits_treemap.png"
    path "RNaquarium_allchunks_blastnanddiamond_allhits_treemap.pdf"
    path "allchunks_alluvialplot_all.png"
    path "allchunks_alluvialplot_all.pdf"

    script:
    """
    mkdir nt
    mkdir nr

    mv $nt_files nt/
    mv $nr_files nr/

    alluvial_treemap.r -nt nt -nr nr
    """
}
