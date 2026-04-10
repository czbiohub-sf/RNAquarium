library(data.table)
library(scales)
library(fs)
library(taxonomizr)
library(tidyverse)
###
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


workingpath <- getwd()
workingpathdash <- str_c(workingpath, "/")
diamondpath <- str_c(workingpath, "/nr_diamond")
blastpath <- str_c(workingpath, "/nt_blast")

## saving all outputs to RNAquarium_outputs
## need to create output folder if it doesn't already exist
dir.create(file.path(workingpath,"RNAquarium_outputs"))

outpath <- str_c(workingpath, "/RNAquarium_outputs")
setwd(outpath)



### combine all NT outputs, then subset viruses
### ADDING AN IF/ELSE SO THAT SCRIPT WILL IMPORT allchunks_blastn_hits_nonhost_ & allchunks_diamond_hits_nonhost_ .tsv.gz files if they already exist!
blastoutputsNTcomplete <- fs::dir_ls(glob="allchunks_blastn_hits*", recurse = FALSE)
blastoutputsNRcomplete <- fs::dir_ls(glob="allchunks_diamond_hits*", recurse = FALSE)



# Process files conditionally
if (length(blastoutputsNTcomplete) == 1 & length(blastoutputsNRcomplete) == 1) {
  # Load existing combined files
  cat("Found existing combined files, loading them...\n")
  thresholded_hit_nofishnomammalsNT <- read_tsv(blastoutputsNTcomplete)
  thresholded_hit_nofishnomammalsNR <- read_tsv(blastoutputsNRcomplete)
  
} else {
  # Combine individual files and save
  cat("Combining individual files...\n")
  
  ## combine all NT outputs
  blastoutputsNT <- fs::dir_ls(blastpath, glob="*.tab", recurse = FALSE)
  thresholded_hit_nofishnomammalsNT <- read_tsv(blastoutputsNT)
  head(thresholded_hit_nofishnomammalsNT)
  write_tsv(thresholded_hit_nofishnomammalsNT, file = paste0("allchunks_blastn_hits_nonhost_",Sys.Date(),".tsv.gz"))
  
  ## combine all NR outputs
  blastoutputsNR <- fs::dir_ls(diamondpath, glob="*.tab", recurse = FALSE)
  thresholded_hit_nofishnomammalsNR <- read_tsv(blastoutputsNR)
  head(thresholded_hit_nofishnomammalsNR)
  write_tsv(thresholded_hit_nofishnomammalsNR, file = paste0("allchunks_diamond_hits_nonhost_",Sys.Date(),".tsv.gz"))
}


########### ADD FULL JOIN TO ALL NOFISHNOMAMMALS, NOT JUST VIRUSES
allchunks_diamondnr_andblastntclustered <- full_join(thresholded_hit_nofishnomammalsNT,thresholded_hit_nofishnomammalsNR)


## first rename!!
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% rename(taxoncategory_NR = taxoncategory2_NR)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% rename(taxoncategory_NTclustered = taxoncategory2_NTclustered)










allchunks_diamondnr_andblastntclustered$taxoncategorysimple_NTclustered <- ifelse((grepl("Root_unresolved", allchunks_diamondnr_andblastntclustered$taxoncategory_NTclustered) == TRUE), "Root_unresolved",
                                                                                  ifelse((grepl("Archaea", allchunks_diamondnr_andblastntclustered$taxoncategory_NTclustered) == TRUE), "Archaea",
                                                                                         ifelse((grepl("Viruses", allchunks_diamondnr_andblastntclustered$taxoncategory_NTclustered) == TRUE), "Viruses",
                                                                                                ifelse((grepl("Bacteria", allchunks_diamondnr_andblastntclustered$taxoncategory_NTclustered) == TRUE), "Bacteria", 
                                                                                                       ifelse((grepl("Fungi", allchunks_diamondnr_andblastntclustered$taxoncategory_NTclustered) == TRUE), "Fungi", 
                                                                                                              ifelse((grepl("Plants", allchunks_diamondnr_andblastntclustered$taxoncategory_NTclustered) == TRUE), "Plants", 
                                                                                                                     ifelse((grepl("SAR_Eukaryotes", allchunks_diamondnr_andblastntclustered$taxoncategory_NTclustered) == TRUE), "SAR_Eukaryotes",
                                                                                                                            ifelse((grepl("other_Eukaryota", allchunks_diamondnr_andblastntclustered$taxoncategory_NTclustered) == TRUE), "other_Eukaryota",
                                                                                                                                   ifelse((grepl("Chordata", allchunks_diamondnr_andblastntclustered$taxoncategory_NTclustered) == TRUE), "Chordata",
                                                                                                                                          ifelse((grepl("Arthropoda", allchunks_diamondnr_andblastntclustered$taxoncategory_NTclustered) == TRUE), "Arthropoda",
                                                                                                                                                 ifelse((grepl("Mollusca", allchunks_diamondnr_andblastntclustered$taxoncategory_NTclustered) == TRUE), "Mollusca",
                                                                                                                                                        ifelse((grepl("Nematoda", allchunks_diamondnr_andblastntclustered$taxoncategory_NTclustered) == TRUE), "Nematoda",
                                                                                                                                                               ifelse((grepl("Annelida", allchunks_diamondnr_andblastntclustered$taxoncategory_NTclustered) == TRUE), "Annelida",
                                                                                                                                                                      ifelse((grepl("Platyhelminthes", allchunks_diamondnr_andblastntclustered$taxoncategory_NTclustered) == TRUE), "Platyhelminthes",
                                                                                                                                                                             ifelse((allchunks_diamondnr_andblastntclustered$taxoncategory_NTclustered == "unresolved_Eukaryota"), "other_Eukaryota", "other_Eukaryota")))))))))))))))

allchunks_diamondnr_andblastntclustered$taxoncategorysimple_NR <- ifelse((grepl("Root_unresolved", allchunks_diamondnr_andblastntclustered$taxoncategory_NR) == TRUE), "Root_unresolved",
                                                                         ifelse((grepl("Archaea", allchunks_diamondnr_andblastntclustered$taxoncategory_NR) == TRUE), "Archaea",
                                                                                ifelse((grepl("Viruses", allchunks_diamondnr_andblastntclustered$taxoncategory_NR) == TRUE), "Viruses",
                                                                                       ifelse((grepl("Bacteria", allchunks_diamondnr_andblastntclustered$taxoncategory_NR) == TRUE), "Bacteria", 
                                                                                              ifelse((grepl("Fungi", allchunks_diamondnr_andblastntclustered$taxoncategory_NR) == TRUE), "Fungi", 
                                                                                                     ifelse((grepl("Plants", allchunks_diamondnr_andblastntclustered$taxoncategory_NR) == TRUE), "Plants", 
                                                                                                            ifelse((grepl("SAR_Eukaryotes", allchunks_diamondnr_andblastntclustered$taxoncategory_NR) == TRUE), "SAR_Eukaryotes",
                                                                                                                   ifelse((grepl("other_Eukaryota", allchunks_diamondnr_andblastntclustered$taxoncategory_NR) == TRUE), "other_Eukaryota",
                                                                                                                          ifelse((grepl("Chordata", allchunks_diamondnr_andblastntclustered$taxoncategory_NR) == TRUE), "Chordata",
                                                                                                                                 ifelse((grepl("Arthropoda", allchunks_diamondnr_andblastntclustered$taxoncategory_NR) == TRUE), "Arthropoda",
                                                                                                                                        ifelse((grepl("Mollusca", allchunks_diamondnr_andblastntclustered$taxoncategory_NR) == TRUE), "Mollusca",
                                                                                                                                               ifelse((grepl("Nematoda", allchunks_diamondnr_andblastntclustered$taxoncategory_NR) == TRUE), "Nematoda",
                                                                                                                                                      ifelse((grepl("Annelida", allchunks_diamondnr_andblastntclustered$taxoncategory_NR) == TRUE), "Annelida",
                                                                                                                                                             ifelse((grepl("Platyhelminthes", allchunks_diamondnr_andblastntclustered$taxoncategory_NR) == TRUE), "Platyhelminthes",
                                                                                                                                                                    ifelse((allchunks_diamondnr_andblastntclustered$taxoncategory_NR == "unresolved_Eukaryota"), "other_Eukaryota", "other_Eukaryota")))))))))))))))



## adding analysis_used, then finally update lowcoverageflag depending on this
allchunks_diamondnr_andblastntclustered <-  allchunks_diamondnr_andblastntclustered %>%
  mutate(
    analysis_used = case_when(
      is.na(bits_NR) ~ "NTclustered",
      is.na(bits_NTclustered) ~ "NR",
      bits_NR >= bits_NTclustered ~ "NR",
      bits_NR < bits_NTclustered ~ "NTclustered",
      .default = "NA"
    )
  )

allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(analysis_used, .before = target_NTclustered)


## finally add the lowcoverage flag only on the analysis_used
allchunks_diamondnr_andblastntclustered <-  allchunks_diamondnr_andblastntclustered %>%
  mutate(
    lowcoverage_flag = case_when(
      (analysis_used == "NR" & sumperc_cov_NR < 0.5) ~ "NR_lowconfidence",
      (analysis_used == "NTclustered" & sumperc_cov_NTclustered < 0.5) ~ "NT_lowconfidence",
      .default = ""
    )
  )

allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(lowcoverage_flag, .before = target_NTclustered)


## MOVE ALL NTORNR FIELDS HERE SO THEY ARE SAVED FOR NEW TAXID_NTORNR SCRIPT
## ALSO FIX MISSPELLING
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, "Japanease", "Japanese")))
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, "Helianthus annus", "Helianthus annuus")))

## MOVE ALL NTORNR FIELDS HERE SO THEY ARE SAVED FOR NEW TAXID_NTORNR SCRIPT
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>%
  mutate(
    taxname_lca_NTorNR = if_else(analysis_used == "NR", taxname_lca_NR, taxname_lca_NTclustered),
    taxoncategory_NTorNR = if_else(analysis_used == "NR", taxoncategory_NR, taxoncategory_NTclustered),
    taxoncategorysimple_NTorNR = if_else(analysis_used == "NR", taxoncategorysimple_NR, taxoncategorysimple_NTclustered),
    bits_NTorNR = if_else(analysis_used == "NR", bits_NR, bits_NTclustered),
    evalue_NTorNR = if_else(analysis_used == "NR", evalue_NR, evalue_NTclustered),
    target_NTorNR = if_else(analysis_used == "NR", target_NR, target_NTclustered),
    taxid_NTorNR = if_else(analysis_used == "NR", taxid_NR, taxid_NTclustered),
    gene_NTorNR = if_else(analysis_used == "NR", gene_NR, gene_NTclustered),
    allele_NTorNR = if_else(analysis_used == "NR", allele_NR, allele_NTclustered),
    pident_NTorNR = if_else(analysis_used == "NR", pident_NR, pident_NTclustered),
    sumperc_cov_NTorNR = if_else(analysis_used == "NR", sumperc_cov_NR, sumperc_cov_NTclustered),
    alnlen_NTorNR = if_else(analysis_used == "NR", alnlen_NR, alnlen_NTclustered),
    mismatch_NTorNR = if_else(analysis_used == "NR", mismatch_NR, mismatch_NTclustered),
    qcov_NTorNR = if_else(analysis_used == "NR", qcov_NR, qcov_NTclustered),
    gapopen_NTorNR = if_else(analysis_used == "NR", gapopen_NR, gapopen_NTclustered),
    qstart_NTorNR = if_else(analysis_used == "NR", qstart_NR, qstart_NTclustered),
    qend_NTorNR = if_else(analysis_used == "NR", qend_NR, qend_NTclustered),
    tstart_NTorNR = if_else(analysis_used == "NR", tstart_NR, tstart_NTclustered),
    tend_NTorNR = if_else(analysis_used == "NR", tend_NR, tend_NTclustered),
    target_title_NTorNR = if_else(analysis_used == "NR", target_title_NR, target_title_NTclustered),
    analysis_NTorNR = if_else(analysis_used == "NR", analysis_NR, analysis_NTclustered),
    sumalnlen_NTorNR = if_else(analysis_used == "NR", sumalnlen_NR, sumalnlen_NTclustered),
    maxbits_NTorNR = if_else(analysis_used == "NR", maxbits_NR, maxbits_NTclustered),
    bits_percmax_NTorNR = if_else(analysis_used == "NR", bits_percmax_NR, bits_percmax_NTclustered),
    tax_superkingdom_NTorNR = if_else(analysis_used == "NR", tax_superkingdom_NR, tax_superkingdom_NTclustered),
    tax_clade_NTorNR = if_else(analysis_used == "NR", tax_clade_NR, tax_clade_NTclustered),
    tax_kingdom_NTorNR = if_else(analysis_used == "NR", tax_kingdom_NR, tax_kingdom_NTclustered),
    tax_phylum_NTorNR = if_else(analysis_used == "NR", tax_phylum_NR, tax_phylum_NTclustered),
    tax_class_NTorNR = if_else(analysis_used == "NR", tax_class_NR, tax_class_NTclustered),
    tax_order_NTorNR = if_else(analysis_used == "NR", tax_order_NR, tax_order_NTclustered),
    tax_family_NTorNR = if_else(analysis_used == "NR", tax_family_NR, tax_family_NTclustered),
    tax_genus_NTorNR = if_else(analysis_used == "NR", tax_genus_NR, tax_genus_NTclustered),
    tax_species_NTorNR = if_else(analysis_used == "NR", tax_species_NR, tax_species_NTclustered)
  )

allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>%
  mutate(
    # First: Replace entries that are just numbers with tax_species_NTorNR
    taxname_lca_NTorNR = ifelse(
      str_detect(taxname_lca_NTorNR, "^[0-9]+$"),  # Pattern for entries that are just numbers
      tax_species_NTorNR,
      taxname_lca_NTorNR
    ),
    # Second: Replace "N/A" entries with tax_species_NTorNR, or tax_genus_NTorNR if species is NA
    taxname_lca_NTorNR = case_when(
      taxname_lca_NTorNR == "N/A" & !is.na(tax_species_NTorNR) ~ tax_species_NTorNR,
      taxname_lca_NTorNR == "N/A" & is.na(tax_species_NTorNR) & !is.na(tax_genus_NTorNR) ~ tax_genus_NTorNR,
      .default = taxname_lca_NTorNR  # Keep original value for all other cases
    ),
    # Third: Replace "N/A" entries with tax_species_NTorNR, or tax_genus_NTorNR if species is NA
    taxname_lca_NTorNR = case_when(
      is.na(taxname_lca_NTorNR) & !is.na(tax_species_NTorNR) ~ tax_species_NTorNR,
      is.na(taxname_lca_NTorNR) & is.na(tax_species_NTorNR) & !is.na(tax_genus_NTorNR) ~ tax_genus_NTorNR,
      .default = taxname_lca_NTorNR  # Keep original value for all other cases
    )
  )

## adding one more catch for even more - but adding a flag!
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>%
  mutate(
    # Fourth: Replace any remaining with taxname_lca_NTclustered or taxname_lca_NR but in these cases add a flag
    lowcoverage_flag = case_when(
      taxname_lca_NTorNR == "N/A" ~ "taxname_lca_lowconfidence",
      is.na(taxname_lca_NTorNR) ~ "taxname_lca_lowconfidence",
      .default = ""
    ),
    taxname_lca_NTorNR = case_when(
      taxname_lca_NTorNR == "N/A" & is.na(taxname_lca_NR) & !is.na(taxname_lca_NTclustered) ~ taxname_lca_NTclustered,
      taxname_lca_NTorNR == "N/A" & taxname_lca_NR == "N/A" & !is.na(taxname_lca_NTclustered) ~ taxname_lca_NTclustered,
      taxname_lca_NTorNR == "N/A" & is.na(taxname_lca_NTclustered) & !is.na(taxname_lca_NR) ~ taxname_lca_NR,
      taxname_lca_NTorNR == "N/A" & taxname_lca_NTclustered == "N/A" & !is.na(taxname_lca_NR) ~ taxname_lca_NR,
      .default = taxname_lca_NTorNR  # Keep original value for all other cases
    ),
    taxname_lca_NTorNR = case_when(
      is.na(taxname_lca_NTorNR) & is.na(taxname_lca_NR) & !is.na(taxname_lca_NTclustered) ~ taxname_lca_NTclustered,
      is.na(taxname_lca_NTorNR) & taxname_lca_NR == "N/A" & !is.na(taxname_lca_NTclustered) ~ taxname_lca_NTclustered,
      is.na(taxname_lca_NTorNR) & is.na(taxname_lca_NTclustered) & !is.na(taxname_lca_NR) ~ taxname_lca_NR,
      is.na(taxname_lca_NTorNR) & taxname_lca_NTclustered == "N/A" & !is.na(taxname_lca_NR) ~ taxname_lca_NR,
      .default = taxname_lca_NTorNR  # Keep original value for all other cases
    )
    
  )

### then relocate
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxname_lca_NTclustered, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategory_NTclustered, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategorysimple_NTclustered, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxname_lca_NR, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategory_NR, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategorysimple_NTclustered, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategorysimple_NR, .before = target_NTclustered)

allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxname_lca_NTorNR, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategory_NTorNR, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategorysimple_NTorNR, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(bits_NTorNR, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(evalue_NTorNR, .before = target_NTclustered)

allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxname_lca_NTclustered, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxname_lca_NR, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategory_NTclustered, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategory_NR, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategorysimple_NTclustered, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategorysimple_NR, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(bits_NTclustered, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(bits_NR, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(evalue_NTclustered, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(evalue_NR, .before = target_NTclustered)

allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(lowcoverage_flag, .before = target_NTclustered)

## adding this to ensure taxname is in 7th column at the end
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(analysis_used, .before = target_NTclustered)


# get mismatches...
allchunks_diamondnr_andblastntclustered_NRbetter <- allchunks_diamondnr_andblastntclustered %>% filter(bits_NR >= bits_NTclustered)
allchunks_diamondnr_andblastntclustered_taxonmismatch <- allchunks_diamondnr_andblastntclustered %>% filter(taxoncategorysimple_NR != taxoncategorysimple_NTclustered)
## then save - NOTE BEFORE REMOVING COLUMNS WE ARE SAVING ONLY THESE 3
## saving full file at at end!!
write_tsv(allchunks_diamondnr_andblastntclustered, file = paste0("taxonomy_hits_nonhost_fullcols_mostrecent.tsv.gz"))

write.table(allchunks_diamondnr_andblastntclustered_NRbetter, file = paste0("taxonomy_hits_nonhost_NRbetter_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_taxonmismatch, file = paste0("taxonomy_hits_nonhost_taxonmismatch_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategory_NTclustered == "Viruses" | taxoncategory_NR == "Viruses")
write.table(allchunks_diamondnr_andblastntclustered_viruses00, file = paste0("taxonomy_hits_viruses0_fullcols_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


## now alluvial plots & treemaps

####################################################################################
###### alluvial code
####################################################################################
allchunks_diamondnr_andblastntclustered_lists <- allchunks_diamondnr_andblastntclustered %>% select(query,taxoncategorysimple_NTorNR,taxoncategorysimple_NTclustered,taxoncategorysimple_NR)

## numbers
## now new alluvial plots...
library(ggalluvial)

data_alluviala <- allchunks_diamondnr_andblastntclustered_lists %>%
  #  group_by(taxoncategorysimple_NTclustered, taxoncategorysimple_NR) %>%
  group_by(taxoncategorysimple_NTclustered, taxoncategorysimple_NTorNR, taxoncategorysimple_NR) %>%
  summarise(count = n()) %>%
  ungroup()


# Replace NA values with 'Missing'
data_alluviala <- data_alluviala %>%
  mutate(across(everything(), ~replace_na(.x, "Missing")))

data_alluviala <- data_alluviala %>% rename(NT = taxoncategorysimple_NTclustered)
data_alluviala <- data_alluviala %>% rename(NR = taxoncategorysimple_NR)
data_alluviala <- data_alluviala %>% rename(used_NTorNR = taxoncategorysimple_NTorNR)

sumcountonepercent <- (sum(data_alluviala$count) / 100)

## set of order of each category using a factor??
data_alluviala <- data_alluviala %>% arrange(desc(count))
alluvialorder2 <- unique(data_alluviala$NT)

## improved version - using all categories but last few gray and not in legend
data_alluvial <- data_alluviala

# Determine category order from ALL data based on used_NTorNR column
# Focus on used_NTorNR since that's what drives the fill aesthetic
used_category_totals <- data_alluvial %>%
  group_by(used_NTorNR) %>%
  summarise(total_count = sum(count), .groups = 'drop') %>%
  arrange(desc(total_count))

# Get all categories present, ordered by abundance, then add "Missing" at the end
all_categories_ordered <- used_category_totals$used_NTorNR
if (!"Missing" %in% all_categories_ordered) {
  all_categories_ordered <- c(all_categories_ordered, "Missing")
}

print(paste("Using", length(all_categories_ordered), "categories in order:"))
print(all_categories_ordered)

# Set up expanded color palette: 11 distinct colors + light gray for the rest
cl_base <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
             "#FFFF33", "#A65628", "#F781BF", "#999999", "#FF62BC", "#A3A500")

# Add light gray colors for categories beyond the first 11
n_categories <- length(all_categories_ordered)
if (n_categories > 11) {
  # Use progressively lighter grays for the remaining categories
  extra_colors <- rep("#D3D3D3", n_categories - 11)  # Light gray for all extras
  cl_expanded <- c(cl_base, extra_colors)
} else {
  cl_expanded <- cl_base[1:n_categories]
}

colors_to_use <- cl_expanded


## for heatmaps below lets use the 11 categories + a new version with all 14
used_category_totalsfull11 <- data_alluvial %>%
  group_by(used_NTorNR) %>%
  summarise(total_count = sum(count), .groups = 'drop') %>%
  arrange(desc(total_count)) %>%
  slice_head(n = 11)

alluvialorder11f <- used_category_totalsfull11$used_NTorNR
alluvialorder11f[12] <- "Missing"

used_category_totalsfull14 <- data_alluvial %>%
  group_by(used_NTorNR) %>%
  summarise(total_count = sum(count), .groups = 'drop') %>%
  arrange(desc(total_count))

alluvialorder14f <- used_category_totalsfull14$used_NTorNR
alluvialorder14f[15] <- "Missing"


## now to figure out 11 hues
cl = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#FF62BC", "#A3A500")
darkcl <- c("#9F1718", "#2D577C", "#397737", "#693870", "#B05601", "#ACAC00", "#743C1E", "#D32092", "#696969", "#CC018C", "#6F7105")
lightcl <- c("#FF6B6B", "#69A4DF", "#6ECD6C", "#C17CCD", "#FFA984", "#FFFF91", "#D4805B", "#FFA7D2", "#B7B7B7", "#FF9ACE", "#BFC144")

## expanded palette with 14 hues
lightcl14 <- c("#FF6B6B", "#69A4DF", "#6ECD6C", "#C17CCD", "#FFA984", "#FFFF91", 
               "#D4805B", "#FFA7D2", "#B7B7B7", "#FF9ACE", "#BFC144",
               "#B3E5FC", "#FFCC80", "#E1BEE7")  # 3 new colors


## edge case when there are groups not in NT, need to combine with other columns
alluvialorder <- data_alluvial %>% pivot_longer(cols = -count)
alluvialorder <- alluvialorder %>% arrange(desc(name),desc(count))
alluvialorder2 <- unique(alluvialorder$value)

# Convert to factors using all categories in abundance order
data_alluvial$NT <- factor(data_alluvial$NT, levels = all_categories_ordered, ordered = TRUE)
data_alluvial$NR <- factor(data_alluvial$NR, levels = all_categories_ordered, ordered = TRUE) 
data_alluvial$used_NTorNR <- factor(data_alluvial$used_NTorNR, levels = all_categories_ordered, ordered = TRUE)


# Create the plot
alluvial_plotall <- ggplot(data_alluvial, aes(axis1 = NT, axis2 = used_NTorNR, axis3 = NR, y = count)) +
  geom_alluvium(aes(fill = used_NTorNR), width = 1/2) +
  stat_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            size = 5, min.y = sumcountonepercent) +
  scale_x_discrete(limits = c("NT", "best (NT or NR)", "NR"), 
                   expand = c(0.08, 0.02)) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = colors_to_use, na.value = "grey50",
                    breaks = all_categories_ordered[1:min(11, length(all_categories_ordered))]) +
  theme_minimal(base_size = 14) +
  labs(title = "Alluvial Diagram of NT vs NR taxonomic categories", 
       fill = "Taxon",
       x = "", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16))

### moving saving of plots to very end

####################################################################################
###### get numbers for treemaps
####################################################################################

### group by then count
data_heatmap <- allchunks_diamondnr_andblastntclustered_lists %>%
  group_by(taxoncategorysimple_NTorNR) %>%
  summarise(count = n()) %>%
  ungroup()

## then sort descending, take just first 10, and use these for treemap
data_heatmap <- data_heatmap %>% arrange(desc(count))
data_heatmap2 <- data_heatmap %>% arrange(desc(count)) %>% slice_head(n = 11)


## new column combining text & numbers, also rewording unite & str_
data_heatmap2 <- data_heatmap2 %>% unite(label, taxoncategorysimple_NTorNR, count, sep = "\n", remove = FALSE)
data_heatmap2 <- data_heatmap2 %>% mutate(label = str_replace_all(label, c("other_Eukaryota" = "other Eukaryota")))
data_heatmap2 <- data_heatmap2 %>% mutate(label = str_replace_all(label, c("Root_unresolved" = "Root (unresolved)")))

library(treemapify)
## 
## updating to get commas in numbers & italics
data_heatmap2 <- data_heatmap2 %>%
  mutate(
    formatted_label = paste0(taxoncategorysimple_NTorNR, "\n", (comma(count)))
  )

## also need to order data_heatmap2
data_heatmap2$taxoncategorysimple_NTorNR <- factor(data_heatmap2$taxoncategorysimple_NTorNR, levels = alluvialorder11f, ordered = TRUE)


NTNRcontigs_treemap <- ggplot2::ggplot(data_heatmap2,aes(area=count,fill=taxoncategorysimple_NTorNR,label=formatted_label,subgroup=taxoncategorysimple_NTorNR)) + 
  treemapify::geom_treemap(layout="squarified") + 
  geom_treemap_text(place = "centre", size = 18, fontface = "italic") + 
  scale_fill_manual(values = lightcl) +
  theme_minimal(base_size = 14) +
  labs(title="Treemap of contigs found by NT + NR searches by taxonomic category", fill="Taxon")


## full 14 version
data_heatmap14 <- data_heatmap %>% arrange(desc(count))

## new column combining text & numbers, also rewording unite & str_
data_heatmap14 <- data_heatmap14 %>% unite(label, taxoncategorysimple_NTorNR, count, sep = "\n", remove = FALSE)
data_heatmap14 <- data_heatmap14 %>% mutate(label = str_replace_all(label, c("other_Eukaryota" = "other Eukaryota")))
data_heatmap14 <- data_heatmap14 %>% mutate(label = str_replace_all(label, c("Root_unresolved" = "Root (unresolved)")))

## updating to get commas in numbers & italics
data_heatmap14 <- data_heatmap14 %>%
  mutate(
    formatted_label = paste0(taxoncategorysimple_NTorNR, "\n", (comma(count)))
  )

## also need to order data_heatmap2
data_heatmap14$taxoncategorysimple_NTorNR <- factor(data_heatmap14$taxoncategorysimple_NTorNR, levels = alluvialorder14f, ordered = TRUE)

NTNRcontigs_treemap14 <- ggplot2::ggplot(data_heatmap14,aes(area=count,fill=taxoncategorysimple_NTorNR,label=formatted_label,subgroup=taxoncategorysimple_NTorNR)) + 
  treemapify::geom_treemap(layout="squarified") + 
  geom_treemap_text(place = "centre", size = 18, fontface = "italic") + 
  scale_fill_manual(values = lightcl14) +
  theme_minimal(base_size = 14) +
  labs(title="Treemap of contigs found by NT + NR searches by taxonomic category", fill="Taxon")


####################################################################################
####################################################################################

## then negative select
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% select(-target_NR) %>% select(-taxid_NTclustered) %>% select(-taxid_NR) %>%
  select(-gene_NTclustered) %>% select(-gene_NR) %>% select(-allele_NTclustered) %>% select(-allele_NR) %>% select(-pident_NTclustered) %>% select(-pident_NR) %>% select(-sumperc_cov_NTclustered) %>% select(-sumperc_cov_NR) %>%
  select(-alnlen_NTclustered) %>% select(-alnlen_NR) %>% select(-mismatch_NTclustered) %>% select(-mismatch_NR) %>% select(-qcov_NTclustered) %>% select(-qcov_NR) %>% select(-gapopen_NTclustered) %>% select(-gapopen_NR) %>%
  select(-qstart_NTclustered) %>% select(-qstart_NR) %>% select(-qend_NTclustered) %>% select(-qend_NR) %>% select(-tstart_NTclustered) %>% select(-tstart_NR) %>% select(-tend_NTclustered) %>% select(-tend_NR) %>% 
  select(-target_title_NTclustered) %>% select(-target_title_NR) %>% select(-analysis_NTclustered) %>% select(-analysis_NR) %>% select(-sumalnlen_NTclustered) %>% select(-sumalnlen_NR) %>% select(-maxbits_NTclustered) %>% select(-maxbits_NR) %>%
  select(-bits_percmax_NTclustered) %>% select(-bits_percmax_NR) %>% select(-tax_superkingdom_NTclustered) %>% select(-tax_superkingdom_NR) %>% select(-tax_clade_NTclustered) %>% select(-tax_clade_NR) %>% select(-tax_kingdom_NTclustered) %>% select(-tax_kingdom_NR) %>%
  select(-tax_phylum_NTclustered) %>% select(-tax_phylum_NR) %>% select(-tax_class_NTclustered) %>% select(-tax_class_NR) %>% select(-tax_order_NTclustered) %>% select(-tax_order_NR) %>% select(-tax_family_NTclustered) %>% select(-tax_family_NR) %>%
  select(-tax_genus_NTclustered) %>% select(-tax_genus_NR) %>% select(-tax_species_NTclustered) %>% select(-tax_species_NR)

## at very end of part 1!
write_tsv(allchunks_diamondnr_andblastntclustered, file = paste0("taxonomy_hits_nonhost_mostrecent.tsv.gz"))


#######################################################################################
#### outputing plots and table last

ggsave(filename = paste("taxonomy_hits_nonhost_alluvialplot_all_",Sys.Date(),".png", sep=""), alluvial_plotall, width = 18, height = 9, units = "in", limitsize = FALSE)
ggsave(filename = paste("taxonomy_hits_nonhost_alluvialplot_all_",Sys.Date(),".pdf", sep=""), alluvial_plotall, width = 18, height = 9, units = "in", limitsize = FALSE)

write.table(data_alluviala, file = paste0("taxonomy_hits_nonhost_alluvialplot_counts_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

ggsave(filename = paste("taxonomy_hits_nonhost_treemap_",Sys.Date(),".png", sep=""), NTNRcontigs_treemap, width = 18, height = 9, units = "in", limitsize = FALSE)
ggsave(filename = paste("taxonomy_hits_nonhost_treemap_",Sys.Date(),".pdf", sep=""), NTNRcontigs_treemap, width = 18, height = 9, units = "in", limitsize = FALSE)

ggsave(filename = paste("taxonomy_hits_nonhost_treemap14categories_",Sys.Date(),".png", sep=""), NTNRcontigs_treemap14, width = 18, height = 9, units = "in", limitsize = FALSE)
ggsave(filename = paste("taxonomy_hits_nonhost_treemap14categories_",Sys.Date(),".pdf", sep=""), NTNRcontigs_treemap14, width = 18, height = 9, units = "in", limitsize = FALSE)


write.table(data_heatmap, file = paste0("taxonomy_hits_nonhost_treemap_counts_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
