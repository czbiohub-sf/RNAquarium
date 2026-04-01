
# Reference index directory names
The default output directory names for mapper indexes are:

 - `hisat2_{Genus_species}_genome`
 - `star_{Genus_species}_indexes.ERCC`
 - `star_{Genus_species}_indexes`
 - `bowtie2_{Genus_species}_index`
 - `gmap_{Genus_species}_genome`
 
(Where `{Genus_species}` comes from the first part of the reference genome filename before a dot
(`.`), as in an ensembl release filename (e.g. `Danio_rerio.GRCz11.dna_sm.primary_assembly.fa`))

These names can be provided to the reference index parameters
 - `--hisat-ref-indexes`
 - `--star-ref-indexes-ercc`
 - `--star-ref-indexes`
 - `--bowtie-ref-indexes`
 - `--gmap-ref-indexes`
respectively, for subsequent runs of the pipeline.
