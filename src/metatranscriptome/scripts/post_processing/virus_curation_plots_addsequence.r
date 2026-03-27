library(data.table)
library(scales)
library(fs)
library(taxonomizr)
library(tidyverse)
library(ggalluvial)
library(phylotools)
library(treemapify)

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
outpath <- str_c(workingpath, "/RNAquarium_outputs")
setwd(outpath)


### read in all files first, then set output path to virus subfolder...
allchunks_diamondnr_andblastntclustered <- read_tsv("taxonomy_hits_viruses0_fullcols_mostrecent.tsv")

## add another if else to check for nonhost_masked.fasta
# Check if masked file exists and load appropriate fasta
if (file.exists("taxonomy_hits_nonhost_list_masked.fasta")) {
  cat("BBDuk masked version of nonhost fasta found, using that fasta for virus sequences\n")
  fasta_viruses <- read.fasta("taxonomy_hits_nonhost_list_masked.fasta")
} else {
  cat("No masked nonhost fasta found, using taxonomy_hits_viruses_list.fasta\n")
  fasta_viruses <- read.fasta("taxonomy_hits_viruses_list.fasta")
}


## need to create output folder if it doesn't already exist
dir.create(file.path(outpath,"virus_outputs"))
outpathvirus <- str_c(outpath, "/virus_outputs")
setwd(outpathvirus)

### AS PART OF VIRUS TAXONOMY - NEED TO FIND TAX_REALM AND REPLACE TAX_CLADE WHICH IS ALL NA FOR VIRUSES...
## this is now placed in step1A - virus curation step that pulls tax_realm using taxonomizr

# Check if tax_clade columns exist and are all NAs - no in viruses00 often occasional Bacteria hits, so needs updating:
if ("tax_clade_NTclustered" %in% names(allchunks_diamondnr_andblastntclustered) && 
    "tax_clade_NR" %in% names(allchunks_diamondnr_andblastntclustered)) {
  
  # Define the expected viral realm names
  viral_realms <- c("Riboviria", "Varidnaviria", "Duplodnaviria")
  
  # Check if any of these realm names exist in either clade column
  realms_in_ntclustered <- any(allchunks_diamondnr_andblastntclustered$tax_clade_NTclustered %in% viral_realms, na.rm = TRUE)
  realms_in_nr <- any(allchunks_diamondnr_andblastntclustered$tax_clade_NR %in% viral_realms, na.rm = TRUE)
  
  if (!realms_in_ntclustered && !realms_in_nr) {
    stop("WARNING: Your virus dataframe is missing the highest level taxonomy, please run script slurm_virus_5A_curation.sh to pull the new 'realm' field to replace your missing values in 'clade'.")
  } else {
    cat("Viral realm taxonomy found - proceeding with analysis.\n")
  }
  
} else {
  # Stop if columns don't exist
  missing_cols <- c("tax_clade_NTclustered", "tax_clade_NR")[!c("tax_clade_NTclustered", "tax_clade_NR") %in% names(allchunks_diamondnr_andblastntclustered)]
  clade_cols <- names(allchunks_diamondnr_andblastntclustered)[grepl("clade", names(allchunks_diamondnr_andblastntclustered), ignore.case = TRUE)]
  
  error_msg <- paste("ERROR: Missing required columns:", paste(missing_cols, collapse = ", "))
  if (length(clade_cols) > 0) {
    error_msg <- paste(error_msg, "\nAvailable clade columns:", paste(clade_cols, collapse = ", "))
  }
  
  stop(error_msg)
}


####################################
## THEN VIRUS SPECIFIC CATEGORIES NEXT

## load 'most recent version' update to use fullcols

allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "Viruses")

## adding very early, replacing all

# Replace "Severe acute respiratory syndrome-related coronavirus" with "Severe acute respiratory syndrome coronavirus 2" across all columns & Betacoronavirus pandemicum
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, "Severe acute respiratory syndrome-related coronavirus", "Severe acute respiratory syndrome coronavirus 2")))
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, "Betacoronavirus pandemicum", "Severe acute respiratory syndrome coronavirus 2")))


## just for viruses, keep a larger set for virus-only alluvial plot
allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategory_NTclustered == "Viruses" | taxoncategory_NR == "Viruses")

# Replace "Severe acute respiratory syndrome-related coronavirus" with "Severe acute respiratory syndrome coronavirus 2" across all columns
allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered_viruses00 %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, "Severe acute respiratory syndrome-related coronavirus", "Severe acute respiratory syndrome coronavirus 2")))
allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered_viruses00 %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, "Betacoronavirus pandemicum", "Severe acute respiratory syndrome coronavirus 2")))

allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, "Japanease", "Japanese")))
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, "Helianthus annus", "Helianthus annuus")))

allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered_viruses00 %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, "Japanease", "Japanese")))
allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered_viruses00 %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, "Helianthus annus", "Helianthus annuus")))

#######################################################################################################################


## now starting to make virus categories

allchunks_diamondnr_andblastntclustered_viruses$viruscategory_NTorNR <- ifelse((grepl("Viruses", allchunks_diamondnr_andblastntclustered_viruses$taxname_lca_NTorNR, ignore.case = TRUE) == TRUE), "Unresolved Viruses",
                                                                               ifelse((grepl("Caudoviricetes", allchunks_diamondnr_andblastntclustered_viruses$tax_class_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                      ifelse((grepl("phage", allchunks_diamondnr_andblastntclustered_viruses$target_title_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                             ifelse((grepl("phix|prokaryotic", allchunks_diamondnr_andblastntclustered_viruses$tax_phylum_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                    ifelse((grepl("phage", allchunks_diamondnr_andblastntclustered_viruses$tax_species_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                           ifelse((grepl("Carnation latent virus|Crocidura shantungensis ribovirus 3|Porcine bastrovirus", allchunks_diamondnr_andblastntclustered_viruses$tax_species_NTclustered, ignore.case = TRUE) == TRUE), "Adapter",
                                                                                                                  ifelse((grepl("Carnation latent virus|Crocidura shantungensis ribovirus 3|Porcine bastrovirus", allchunks_diamondnr_andblastntclustered_viruses$taxname_lca_NTclustered, ignore.case = TRUE) == TRUE), "Adapter",
                                                                                                                         ifelse((grepl("Porcine picobirnavirus", allchunks_diamondnr_andblastntclustered_viruses$tax_species_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                ifelse((grepl("Levivir", allchunks_diamondnr_andblastntclustered_viruses$tax_class_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                       ifelse((grepl("Inoviridae|Levivir", allchunks_diamondnr_andblastntclustered_viruses$tax_family_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                              ifelse((grepl("Cystoviridae", allchunks_diamondnr_andblastntclustered_viruses$tax_family_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                                     ifelse((grepl("phage|Prokaryotic", allchunks_diamondnr_andblastntclustered_viruses$taxname_lca_NTorNR, ignore.case = TRUE) == TRUE), "Phage", allchunks_diamondnr_andblastntclustered_viruses$taxname_lca_NTorNR))))))))))))

allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(viruscategory_NTorNR, .before = taxname_lca_NTclustered)

allchunks_diamondnr_andblastntclustered_viruses00$viruscategory_NTorNR <- ifelse((grepl("Viruses", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NTorNR, ignore.case = TRUE) == TRUE), "Unresolved Viruses",
                                                                                 ifelse((grepl("Caudoviricetes", allchunks_diamondnr_andblastntclustered_viruses00$tax_class_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                        ifelse((grepl("phage", allchunks_diamondnr_andblastntclustered_viruses00$target_title_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                               ifelse((grepl("phix|prokaryotic", allchunks_diamondnr_andblastntclustered_viruses00$tax_phylum_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                      ifelse((grepl("phage", allchunks_diamondnr_andblastntclustered_viruses00$tax_species_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                             ifelse((grepl("Carnation latent virus|Crocidura shantungensis ribovirus 3|Porcine bastrovirus", allchunks_diamondnr_andblastntclustered_viruses00$tax_species_NTclustered, ignore.case = TRUE) == TRUE), "Adapter",
                                                                                                                    ifelse((grepl("Carnation latent virus|Crocidura shantungensis ribovirus 3|Porcine bastrovirus", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NTclustered, ignore.case = TRUE) == TRUE), "Adapter",
                                                                                                                           ifelse((grepl("Porcine picobirnavirus", allchunks_diamondnr_andblastntclustered_viruses00$tax_species_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                  ifelse((grepl("Levivir", allchunks_diamondnr_andblastntclustered_viruses$tax_class_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                         ifelse((grepl("Inoviridae|Levivir", allchunks_diamondnr_andblastntclustered_viruses$tax_family_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                                ifelse((grepl("Cystoviridae", allchunks_diamondnr_andblastntclustered_viruses$tax_family_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                                       ifelse((grepl("phage|Prokaryotic", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NTorNR, ignore.case = TRUE) == TRUE), "Phage", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NTorNR))))))))))))

allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered_viruses00 %>% relocate(viruscategory_NTorNR, .before = taxname_lca_NTclustered)




allchunks_diamondnr_andblastntclustered_viruses$viruscategorysimple_NTorNR <- ifelse((grepl("Viruses", allchunks_diamondnr_andblastntclustered_viruses$taxname_lca_NTorNR, ignore.case = TRUE) == TRUE), "Unresolved Viruses",
                                                                                     ifelse((grepl("Caudoviricetes", allchunks_diamondnr_andblastntclustered_viruses$tax_class_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                            ifelse((grepl("phage", allchunks_diamondnr_andblastntclustered_viruses$target_title_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                   ifelse((grepl("phix|prokaryotic", allchunks_diamondnr_andblastntclustered_viruses$tax_phylum_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                          ifelse((grepl("phage", allchunks_diamondnr_andblastntclustered_viruses$tax_species_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                 ifelse((grepl("Carnation latent virus|Crocidura shantungensis ribovirus 3|Porcine bastrovirus", allchunks_diamondnr_andblastntclustered_viruses$tax_species_NTclustered, ignore.case = TRUE) == TRUE), "Adapter",
                                                                                                                        ifelse((grepl("Carnation latent virus|Crocidura shantungensis ribovirus 3|Porcine bastrovirus", allchunks_diamondnr_andblastntclustered_viruses$taxname_lca_NTclustered, ignore.case = TRUE) == TRUE), "Adapter",
                                                                                                                               ifelse((grepl("Porcine picobirnavirus", allchunks_diamondnr_andblastntclustered_viruses$tax_species_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                      ifelse((grepl("Levivir", allchunks_diamondnr_andblastntclustered_viruses$tax_class_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                             ifelse((grepl("Inoviridae|Levivir", allchunks_diamondnr_andblastntclustered_viruses$tax_family_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                                    ifelse((grepl("Cystoviridae", allchunks_diamondnr_andblastntclustered_viruses$tax_family_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                                           ifelse((grepl("phage|Prokaryotic", allchunks_diamondnr_andblastntclustered_viruses$taxname_lca_NTorNR, ignore.case = TRUE) == TRUE), "Phage", "Non-phage"))))))))))))

allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(viruscategorysimple_NTorNR, .before = taxname_lca_NTclustered)


## taxoncategorysimple_NTorNR =! "Viruses"
allchunks_diamondnr_andblastntclustered_viruses00$viruscategorysimple_NTorNR <- ifelse((grepl("Viruses", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NTorNR, ignore.case = TRUE) == TRUE), "Unresolved Viruses",
                                                                                       ifelse((grepl("Caudoviricetes", allchunks_diamondnr_andblastntclustered_viruses00$tax_class_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                              ifelse((grepl("phage", allchunks_diamondnr_andblastntclustered_viruses00$target_title_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                     ifelse((grepl("phix|prokaryotic", allchunks_diamondnr_andblastntclustered_viruses00$tax_phylum_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                            ifelse((grepl("phage", allchunks_diamondnr_andblastntclustered_viruses00$tax_species_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                   ifelse((grepl("Carnation latent virus|Crocidura shantungensis ribovirus 3|Porcine bastrovirus", allchunks_diamondnr_andblastntclustered_viruses00$tax_species_NTclustered, ignore.case = TRUE) == TRUE), "Adapter",
                                                                                                                          ifelse((grepl("Carnation latent virus|Crocidura shantungensis ribovirus 3|Porcine bastrovirus", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NTclustered, ignore.case = TRUE) == TRUE), "Adapter",
                                                                                                                                 ifelse((grepl("Porcine picobirnavirus", allchunks_diamondnr_andblastntclustered_viruses00$tax_species_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                        ifelse((grepl("Levivir", allchunks_diamondnr_andblastntclustered_viruses$tax_class_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                               ifelse((grepl("Inoviridae|Levivir", allchunks_diamondnr_andblastntclustered_viruses$tax_family_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                                      ifelse((grepl("Cystoviridae", allchunks_diamondnr_andblastntclustered_viruses$tax_family_NTorNR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                                             ifelse((grepl("phage|Prokaryotic", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NTorNR, ignore.case = TRUE) == TRUE), "Phage", "Non-phage"))))))))))))

allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered_viruses00 %>% relocate(viruscategorysimple_NTorNR, .before = taxname_lca_NTclustered)

allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered_viruses00 %>%
  mutate(
    viruscategorysimple_NTorNR = if_else(taxoncategorysimple_NTorNR != "Viruses", "Non-virus", viruscategorysimple_NTorNR)
  )

## we actually need the old columns (just NT & just NR) for alluvial plots in viruses00 file only
## careful - do not use non-phage as default if there are missing values, also adapters are too few
## moving missing up to here
allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered_viruses00 %>%
  mutate(across(c(taxname_lca_NR, taxname_lca_NTclustered), ~replace_na(.x, "Missing")))

allchunks_diamondnr_andblastntclustered_viruses01 <- allchunks_diamondnr_andblastntclustered_viruses00


allchunks_diamondnr_andblastntclustered_viruses00$viruscategorysimple_NR <- ifelse((grepl("Viruses", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NR, ignore.case = TRUE) == TRUE), "Unresolved Viruses",
                                                                                   ifelse((grepl("Caudoviricetes", allchunks_diamondnr_andblastntclustered_viruses00$tax_class_NR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                          ifelse((grepl("phage", allchunks_diamondnr_andblastntclustered_viruses00$target_title_NR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                 ifelse((grepl("phix|prokaryotic", allchunks_diamondnr_andblastntclustered_viruses00$tax_phylum_NR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                        ifelse((grepl("phage", allchunks_diamondnr_andblastntclustered_viruses00$tax_species_NR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                               ifelse((grepl("Crocidura shantungensis ribovirus 3|Porcine bastrovirus", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NR, ignore.case = TRUE) == TRUE), "Adapter",
                                                                                                                      ifelse((grepl("Levivir", allchunks_diamondnr_andblastntclustered_viruses00$tax_class_NR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                             ifelse((grepl("Inoviridae", allchunks_diamondnr_andblastntclustered_viruses00$tax_family_NR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                    ifelse((grepl("Cystoviridae|Levivir", allchunks_diamondnr_andblastntclustered_viruses00$tax_family_NR, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                           ifelse((grepl("Missing", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NR, ignore.case = TRUE) == TRUE), "Missing",
                                                                                                                                                  ifelse((grepl("phage|Prokaryotic|Porcine picobirnavirus", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NR, ignore.case = TRUE) == TRUE), "Phage", "Non-phage")))))))))))


allchunks_diamondnr_andblastntclustered_viruses00$viruscategorysimple_NTclustered <- ifelse((grepl("Viruses", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NTclustered, ignore.case = TRUE) == TRUE), "Unresolved Viruses",
                                                                                            ifelse((grepl("Caudoviricetes", allchunks_diamondnr_andblastntclustered_viruses00$tax_class_NTclustered, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                   ifelse((grepl("phage", allchunks_diamondnr_andblastntclustered_viruses00$target_title_NTclustered, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                          ifelse((grepl("phix|prokaryotic", allchunks_diamondnr_andblastntclustered_viruses00$tax_phylum_NTclustered, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                 ifelse((grepl("phage", allchunks_diamondnr_andblastntclustered_viruses00$tax_species_NTclustered, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                        ifelse((grepl("Carnation latent virus|Crocidura shantungensis ribovirus 3|Porcine bastrovirus", allchunks_diamondnr_andblastntclustered_viruses00$tax_species_NTclustered, ignore.case = TRUE) == TRUE), "Adapter",
                                                                                                                               ifelse((grepl("Carnation latent virus|Crocidura shantungensis ribovirus 3|Porcine bastrovirus", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NTclustered, ignore.case = TRUE) == TRUE), "Adapter",
                                                                                                                                      ifelse((grepl("Porcine picobirnavirus", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NTclustered, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                             ifelse((grepl("Levivir", allchunks_diamondnr_andblastntclustered_viruses00$tax_class_NTclustered, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                                    ifelse((grepl("Inoviridae|Levivir", allchunks_diamondnr_andblastntclustered_viruses00$tax_family_NTclustered, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                                           ifelse((grepl("Cystoviridae", allchunks_diamondnr_andblastntclustered_viruses00$tax_family_NTclustered, ignore.case = TRUE) == TRUE), "Phage",
                                                                                                                                                                  ifelse((grepl("Missing", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NTclustered, ignore.case = TRUE) == TRUE), "Missing",
                                                                                                                                                                         ifelse((grepl("phage|Prokaryotic", allchunks_diamondnr_andblastntclustered_viruses00$taxname_lca_NTclustered, ignore.case = TRUE) == TRUE), "Phage", "Non-phage")))))))))))))




allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered_viruses00 %>% relocate(viruscategorysimple_NR, .before = taxname_lca_NTclustered)
allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered_viruses00 %>% relocate(viruscategorysimple_NTclustered, .before = taxname_lca_NTclustered)


## first need this:
allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered_viruses00 %>%
  mutate(across(c(taxoncategory_NR, taxoncategory_NTclustered), ~replace_na(.x, "Missing")))
allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered_viruses00 %>%
  mutate(across(c(taxoncategorysimple_NR, taxoncategorysimple_NTclustered), ~replace_na(.x, "Missing")))

allchunks_diamondnr_andblastntclustered_viruses %>% group_by(viruscategorysimple_NTorNR) %>% summarize(count=n())


## edge case where bioproject isnt named
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>%
  mutate(bioproject = if_else(bioproject == "NA", "other", bioproject))
allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered_viruses00 %>%
  mutate(bioproject = if_else(bioproject == "NA", "other", bioproject))


## then save only files before removing columns here, then rest at end

write.table(allchunks_diamondnr_andblastntclustered_viruses, file = paste0("taxonomy_hits_viruses_fullcols_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


## also save mismatches for viruses
# get mismatches...
allchunks_diamondnr_andblastntclustered_viruses00_NRbetter <- allchunks_diamondnr_andblastntclustered_viruses00 %>% filter(bits_NR >= bits_NTclustered)
allchunks_diamondnr_andblastntclustered_viruses00_taxonmismatch <- allchunks_diamondnr_andblastntclustered_viruses00 %>% filter(taxoncategorysimple_NR != taxoncategorysimple_NTclustered)
## then save
write.table(allchunks_diamondnr_andblastntclustered_viruses00_NRbetter, file = paste0("taxonomy_hits_viruses00_NRbetter_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_viruses00_taxonmismatch, file = paste0("taxonomy_hits_viruses00_taxonmismatch_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


###############################################################
####### NOW ALLUVIAL PLOTS - AGAIN HERE USING VIRUSES00 VERSION...


## UPDATE - INCLUDE MIDDLE COLUMN OF ACTUALLY USED...ALSO REMOVE NON-VIRUSES BEFORE ALLUVIAL
data_alluvialv <- allchunks_diamondnr_andblastntclustered_viruses00 %>%
  group_by(viruscategorysimple_NTclustered, viruscategorysimple_NTorNR, viruscategorysimple_NR) %>%
  summarise(count = n()) %>%
  ungroup()

data_alluvialv <- data_alluvialv %>% rename(NT = viruscategorysimple_NTclustered)
data_alluvialv <- data_alluvialv %>% rename(NR = viruscategorysimple_NR)
data_alluvialv <- data_alluvialv %>% rename(used_NTorNR = viruscategorysimple_NTorNR)

## filter non-virus
data_alluvialv <- data_alluvialv %>% dplyr::filter(used_NTorNR != "Non-virus")
sumcountonepercent <- (sum(data_alluvialv$count) / 100)


## set of order of each category using a factor??
data_alluvialv <- data_alluvialv %>% arrange(desc(count))
alluvialorder2 <- unique(data_alluvialv$NT)
## edge case when there are groups not in NT, need to combine with other columns
alluvialorder <- data_alluvialv %>% pivot_longer(cols = -count)
alluvialorder <- alluvialorder %>% arrange(desc(name),desc(count))
alluvialorder2 <- unique(alluvialorder$value)


data_alluvialv$NT <- factor(data_alluvialv$NT, levels = alluvialorder2, ordered = TRUE)
data_alluvialv$NR <- factor(data_alluvialv$NR, levels = alluvialorder2, ordered = TRUE)
data_alluvialv$used_NTorNR <- factor(data_alluvialv$used_NTorNR, levels = alluvialorder2, ordered = TRUE)
## this works, but also have to remove "decreasing = FALSE" across the commands

## we will want to use alluvialorder2 order also in treemap down below
if ("Adapter" %in% unique(data_alluvialv$used_NTorNR)) {
  clv <- c("#00969A", "#4D6E00", "#C12600", "#9E00DF")
} else {
  # color_values <- c("#7CAE00", "#00BFC4", "#C77CFF")
  clv <- c("#00969A", "#4D6E00", "#9E00DF")
}


alluvial_plotv <- ggplot(data_alluvialv, aes(axis1 = NT, axis2 = used_NTorNR, axis3 = NR, y = count)) +
  geom_alluvium(aes(fill = used_NTorNR), width = 1/2) +
  stat_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 5, min.y = sumcountonepercent) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(limits = c("NT", "best (NT or NR)", "NR"), expand = c(0.08, 0.02)) +
  theme_minimal(base_size = 14) + ## also theme_minimal, theme_void
  #  theme_minimal(base_family="Helvetica", base_size = 16) + ## also theme_minimal, theme_void
  # scale_fill_brewer(palette = "Dark2") + 
  # scale_fill_manual(values = c("#00BFC4", "#7CAE00" ,"#F8766D","#C77CFF")) + 
  #  scale_fill_manual(values = c("#058286", "#547700", "#CE2C0A", "#AA00EF")) + ## darkened (see below) by 0.3
  #  scale_fill_manual(values = c("#00969A", "#618907", "#CF544B", "#A854E0")) + ## 0.2
  #  scale_fill_manual(values = c( "#00969A", "#4D6E00", "#C12600", "#9E00DF")) + ## 0.35 can use names like values = darkcl
  scale_fill_manual(values = clv) + ## 0.35 can use names like values = darkcl
  labs(title = "Alluvial Diagram of NT vs NR virus categories", fill = "Virus Category",
       x = "",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16))

############################################################################
## now remove unused columns...no need to do this for viruses00 file, since that isn't saved
## do just for viruses then create interesting set after
## then also consolidate many columns, then just remove a bunch

## we added a cleaner version above and on more columns!

## then negative select
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% select(-target_NR) %>% select(-taxid_NTclustered) %>% select(-taxid_NR) %>%
  select(-gene_NTclustered) %>% select(-gene_NR) %>% select(-allele_NTclustered) %>% select(-allele_NR) %>% select(-pident_NTclustered) %>% select(-pident_NR) %>% select(-sumperc_cov_NTclustered) %>% select(-sumperc_cov_NR) %>%
  select(-alnlen_NTclustered) %>% select(-alnlen_NR) %>% select(-mismatch_NTclustered) %>% select(-mismatch_NR) %>% select(-qcov_NTclustered) %>% select(-qcov_NR) %>% select(-gapopen_NTclustered) %>% select(-gapopen_NR) %>%
  select(-qstart_NTclustered) %>% select(-qstart_NR) %>% select(-qend_NTclustered) %>% select(-qend_NR) %>% select(-tstart_NTclustered) %>% select(-tstart_NR) %>% select(-tend_NTclustered) %>% select(-tend_NR) %>% 
  select(-target_title_NTclustered) %>% select(-target_title_NR) %>% select(-analysis_NTclustered) %>% select(-analysis_NR) %>% select(-sumalnlen_NTclustered) %>% select(-sumalnlen_NR) %>% select(-maxbits_NTclustered) %>% select(-maxbits_NR) %>%
  select(-bits_percmax_NTclustered) %>% select(-bits_percmax_NR) %>% select(-tax_superkingdom_NTclustered) %>% select(-tax_superkingdom_NR) %>% select(-tax_clade_NTclustered) %>% select(-tax_clade_NR) %>% select(-tax_kingdom_NTclustered) %>% select(-tax_kingdom_NR) %>%
  select(-tax_phylum_NTclustered) %>% select(-tax_phylum_NR) %>% select(-tax_class_NTclustered) %>% select(-tax_class_NR) %>% select(-tax_order_NTclustered) %>% select(-tax_order_NR) %>% select(-tax_family_NTclustered) %>% select(-tax_family_NR) %>%
  select(-tax_genus_NTclustered) %>% select(-tax_genus_NR) %>% select(-tax_species_NTclustered) %>% select(-tax_species_NR)



## also for viruses change some taxonommy to NAs when seeing "viria" "virae" "viricota" "viricetes" "virales" "viridae"
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>%
  mutate(
    # For rows where 'taxname_lca_NTorNR' is just Viruses
    tax_clade_NTorNR = if_else(grepl("Viruses$", taxname_lca_NTorNR), NA, tax_clade_NTorNR),
    tax_kingdom_NTorNR = if_else(grepl("Viruses$", taxname_lca_NTorNR), NA, tax_kingdom_NTorNR),
    tax_phylum_NTorNR = if_else(grepl("Viruses$", taxname_lca_NTorNR), NA, tax_phylum_NTorNR),
    tax_class_NTorNR = if_else(grepl("Viruses$", taxname_lca_NTorNR), NA, tax_class_NTorNR),
    tax_order_NTorNR = if_else(grepl("Viruses$", taxname_lca_NTorNR), NA, tax_order_NTorNR),
    tax_family_NTorNR = if_else(grepl("Viruses$", taxname_lca_NTorNR), NA, tax_family_NTorNR),
    tax_genus_NTorNR = if_else(grepl("Viruses$", taxname_lca_NTorNR), NA, tax_genus_NTorNR),
    
    # For rows where 'taxname_lca_NTorNR' ends with "viria"
    tax_kingdom_NTorNR = if_else(grepl("viria$", taxname_lca_NTorNR), NA, tax_kingdom_NTorNR),
    tax_phylum_NTorNR = if_else(grepl("viria$", taxname_lca_NTorNR), NA, tax_phylum_NTorNR),
    tax_class_NTorNR = if_else(grepl("viria$", taxname_lca_NTorNR), NA, tax_class_NTorNR),
    tax_order_NTorNR = if_else(grepl("viria$", taxname_lca_NTorNR), NA, tax_order_NTorNR),
    tax_family_NTorNR = if_else(grepl("viria$", taxname_lca_NTorNR), NA, tax_family_NTorNR),
    tax_genus_NTorNR = if_else(grepl("viria$", taxname_lca_NTorNR), NA, tax_genus_NTorNR),
    
    # For rows where 'taxname_lca_NTorNR' ends with "virae"
    tax_phylum_NTorNR = if_else(grepl("virae$", taxname_lca_NTorNR), NA, tax_phylum_NTorNR),
    tax_class_NTorNR = if_else(grepl("virae$", taxname_lca_NTorNR), NA, tax_class_NTorNR),
    tax_order_NTorNR = if_else(grepl("virae$", taxname_lca_NTorNR), NA, tax_order_NTorNR),
    tax_family_NTorNR = if_else(grepl("virae$", taxname_lca_NTorNR), NA, tax_family_NTorNR),
    tax_genus_NTorNR = if_else(grepl("virae$", taxname_lca_NTorNR), NA, tax_genus_NTorNR),
    
    # For rows where 'taxname_lca_NTorNR' ends with "viricota"
    tax_class_NTorNR = if_else(grepl("viricota$", taxname_lca_NTorNR), NA, tax_class_NTorNR),
    tax_order_NTorNR = if_else(grepl("viricota$", taxname_lca_NTorNR), NA, tax_order_NTorNR),
    tax_family_NTorNR = if_else(grepl("viricota$", taxname_lca_NTorNR), NA, tax_family_NTorNR),
    tax_genus_NTorNR = if_else(grepl("viricota$", taxname_lca_NTorNR), NA, tax_genus_NTorNR),
    
    # For rows where 'taxname_lca_NTorNR' ends with "viricetes"
    tax_order_NTorNR = if_else(grepl("viricetes$", taxname_lca_NTorNR), NA, tax_order_NTorNR),
    tax_family_NTorNR = if_else(grepl("viricetes$", taxname_lca_NTorNR), NA, tax_family_NTorNR),
    tax_genus_NTorNR = if_else(grepl("viricetes$", taxname_lca_NTorNR), NA, tax_genus_NTorNR),
    
    # For rows where 'taxname_lca_NTorNR' ends with "virales"
    tax_family_NTorNR = if_else(grepl("virales$", taxname_lca_NTorNR), NA, tax_family_NTorNR),
    tax_genus_NTorNR = if_else(grepl("virales$", taxname_lca_NTorNR), NA, tax_genus_NTorNR),
    
    # For rows where 'taxname_lca_NTorNR' ends with "viridae"
    tax_genus_NTorNR = if_else(grepl("viridae$", taxname_lca_NTorNR), NA, tax_genus_NTorNR)
  )

## relocate taxonomy columns
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(tax_species_NTorNR, .after = target_title_NTorNR)
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(tax_genus_NTorNR, .after = target_title_NTorNR)
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(tax_family_NTorNR, .after = target_title_NTorNR)
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(tax_order_NTorNR, .after = target_title_NTorNR)
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(tax_class_NTorNR, .after = target_title_NTorNR)
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(tax_phylum_NTorNR, .after = target_title_NTorNR)
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(tax_kingdom_NTorNR, .after = target_title_NTorNR)
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(tax_clade_NTorNR, .after = target_title_NTorNR)
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(tax_superkingdom_NTorNR, .after = target_title_NTorNR)


#allchunks_diamondnr_andblastntclustered_viruses <- read_tsv("taxonomy_hits_viruses_mostrecent.tsv")

## more relocates
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(taxname_lca_NTclustered, .after = last_col())
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(taxname_lca_NR, .after = last_col())
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(target_NTclustered, .after = last_col())
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(taxoncategory_NTclustered, .after = last_col())
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(taxoncategory_NR, .after = last_col())
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(taxoncategorysimple_NTclustered, .after = last_col())
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(taxoncategorysimple_NR, .after = last_col())
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(bits_NTclustered, .after = last_col())
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(bits_NR, .after = last_col())
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(evalue_NTclustered, .after = last_col())
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(evalue_NR, .after = last_col())

allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(lowcoverage_flag, .before = taxid_NTorNR)
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(target_title_NTorNR, .before = taxid_NTorNR)
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(target_NTorNR, .before = taxid_NTorNR)



## then finally create the more interesting subset - filter out phage, adapter 
allchunks_diamondnr_andblastntclustered_nonphage <- allchunks_diamondnr_andblastntclustered_viruses
allchunks_diamondnr_andblastntclustered_nonphage <- allchunks_diamondnr_andblastntclustered_nonphage %>% dplyr::filter(viruscategorysimple_NTorNR != "Unresolved Viruses")
allchunks_diamondnr_andblastntclustered_nonphage <- allchunks_diamondnr_andblastntclustered_nonphage %>% dplyr::filter(viruscategorysimple_NTorNR != "Adapter")
allchunks_diamondnr_andblastntclustered_nonphage <- allchunks_diamondnr_andblastntclustered_nonphage %>% dplyr::filter(viruscategorysimple_NTorNR != "Phage")


## now save ONLY viruses & interesting versions - viruses00 is only used for alluvial plots
write.table(allchunks_diamondnr_andblastntclustered_viruses, file = paste0("taxonomy_hits_viruses_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

write.table(allchunks_diamondnr_andblastntclustered_nonphage, file = paste0("taxonomy_hits_viruses_nonphage_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## also always have a generic version saved NOTE THIS WILL CHANGE, BUT IS USED FOR PART 3 AND ALSO VIRUS CURATION...
write.table(allchunks_diamondnr_andblastntclustered_nonphage, file = paste0("taxonomy_hits_viruses_nonphage_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

write.table(allchunks_diamondnr_andblastntclustered_viruses, file = paste0("taxonomy_hits_viruses_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

#############################################################################
### ALSO ADD SEQUENCE JUST AS IN PART3, BUT FOR NEWLY CREATED VIRUS FILES...

## rename seq.name to fullquery
fasta_viruses <- fasta_viruses %>% rename(query = seq.name)
allchunks_diamondnr_andblastntclustered_viruses_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_viruses, fasta_viruses)
write.table(allchunks_diamondnr_andblastntclustered_viruses_withsequence, file = paste0("taxonomy_hits_hits_viruses_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_nonphage, fasta_viruses)
write.table(allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence, file = paste0("taxonomy_hits_hits_viruses_nonphage_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## again always have a generic version saved for after clustering
write.table(allchunks_diamondnr_andblastntclustered_viruses_withsequence, file = paste0("taxonomy_hits_viruses_withsequence_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence, file = paste0("taxonomy_hits_viruses_nonphage_withsequence_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


## let's also take the most critical columns, add their info to the query and output fastas here
## columns to include

## use unite() on many columns using _ as sep, then calling it 'read', then use writetoFasta function...
## columns to include - query, taxname_lca_NTorNR, analysis_used, evalue_NTorNR, lowcoverage_flag
allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence %>% unite(read, query, taxname_lca_NTorNR, analysis_used, evalue_NTorNR, lowcoverage_flag, sep = "|", remove = FALSE)
allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence_forfasta %>% select(read,seq.text)
## drop "_NA" with gsub
allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence_forfasta$read <- gsub('_lowconfidence','lowconfidence',allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence_forfasta$read)
allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence_forfasta$read <- gsub(' ','_',allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence_forfasta$read)


## NEED TO REMOVE SPACES...HOW ABOUT FIRST CHANGE _ TO | THEN SPACES TO _
allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence_forfasta$read <- gsub(' ','_',allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence_forfasta$read)


writetoFasta <- function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"read"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq.text"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}
## optimized
writetoFastafaster <- function(data, filename) {
  # Construct the FASTA lines in a vectorized way
  fastaLines <- paste0(">", data[["read"]], "\n", data[["seq.text"]])
  
  # Write directly to the file
  writeLines(fastaLines, filename)
}


writetoFastafaster(allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence_forfasta, paste0("taxonomy_hits_viruses_nonphage_",Sys.Date(),".fasta"))


## also run for all viruses, not just interesting - note old fasta list allchunks_blastnanddiamond_hits_viruses_list.fasta, new is allchunks_blastnanddiamond_hits_viruses.fasta & allchunks_blastnanddiamond_hits_viruses_nonphage.fasta
allchunks_diamondnr_andblastntclustered_viruses_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_withsequence %>% unite(read, query, taxname_lca_NTorNR, analysis_used, evalue_NTorNR, lowcoverage_flag, sep = "|", remove = FALSE)
allchunks_diamondnr_andblastntclustered_viruses_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_withsequence_forfasta %>% select(read,seq.text)
## drop "_NA" with gsub
allchunks_diamondnr_andblastntclustered_viruses_withsequence_forfasta$read <- gsub('_lowconfidence','lowconfidence',allchunks_diamondnr_andblastntclustered_viruses_withsequence_forfasta$read)

allchunks_diamondnr_andblastntclustered_viruses_withsequence_forfasta$read <- gsub(' ','_',allchunks_diamondnr_andblastntclustered_viruses_withsequence_forfasta$read)


writetoFastafaster(allchunks_diamondnr_andblastntclustered_viruses_withsequence_forfasta, paste0("taxonomy_hits_viruses_",Sys.Date(),".fasta"))

## also mostrecent not just with dates
writetoFastafaster(allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence_forfasta, "taxonomy_hits_viruses_nonphage_mostrecent.fasta")
writetoFastafaster(allchunks_diamondnr_andblastntclustered_viruses_withsequence_forfasta, "taxonomy_hits_viruses_mostrecent.fasta")

####################################################################################
###### get numbers for treemaps
####################################################################################

## UPDATE - INCLUDE COUNT OF BIOPROJECTS BY LCA
virus_heatmap0 <- allchunks_diamondnr_andblastntclustered_nonphage %>% group_by(taxname_lca_NTorNR) %>% mutate(bioprojectcount = n_distinct(bioproject, na.rm = TRUE)) %>% ungroup()

virus_heatmap <- virus_heatmap0 %>%
  group_by(taxname_lca_NTorNR) %>%
  summarise(count = n(), bioprojectcount = first(bioprojectcount)) %>%
  ungroup()


## then sort descending, take just first 25, and use these for treemap
virus_heatmap <- virus_heatmap %>% arrange(desc(count),desc(bioprojectcount))
virus_heatmap2 <- virus_heatmap %>% arrange(desc(count),desc(bioprojectcount)) %>% slice_head(n = 25)

## save these numbers (moving to end)

## new column combining text & numbers, also rewording unite & str_ Severe acute respiratory syndrome coronavirus 2
virus_heatmap2 <- virus_heatmap2 %>% mutate(taxname_lca_NTorNR = str_replace_all(taxname_lca_NTorNR, c("Severe acute respiratory syndrome coronavirus 2" = "SARS-CoV-2")))
virus_heatmap2 <- virus_heatmap2 %>% mutate(taxname_lca_NTorNR = str_replace_all(taxname_lca_NTorNR, c("Human immunodeficiency virus 1" = "HIV-1")))

## changing tab below to space - tab breaks in png version...
virus_heatmap2 <- virus_heatmap2 %>%
  mutate(
    formatted_label = paste0(taxname_lca_NTorNR, "\n", (comma(count)), " ", (comma(bioprojectcount)))
  )


virus_treemap <- ggplot2::ggplot(virus_heatmap2,aes(area=count,fill=taxname_lca_NTorNR,label=formatted_label,subgroup=taxname_lca_NTorNR)) + 
  treemapify::geom_treemap(layout="squarified") + 
  geom_treemap_text(place = "centre", size = 18, fontface = "italic") + 
  theme_minimal(base_size = 14, base_family = "sans") +
  labs(title="Treemap of non-phage viral contigs found by NT + NR searches, with contig & bioproject counts per LCA", fill="Taxon") +
  theme(legend.position = "none")  # Removes the external legend

## saving all plots at very end

###################################
### another treemap showing broad virus categories:
virus_heatmapbroad <- allchunks_diamondnr_andblastntclustered_viruses %>% mutate(subgroup = ifelse(viruscategorysimple_NTorNR == "Non-phage", tax_clade_NTorNR, viruscategorysimple_NTorNR))

#virus_heatmapbroad <- virus_heatmapbroad %>% group_by(viruscategorysimple_NTorNR,subgroup)
virus_heatmapbroad <- virus_heatmapbroad %>%
  group_by(viruscategorysimple_NTorNR,subgroup) %>%
  summarise(count = n()) %>%
  ungroup()

virus_heatmapbroad <- virus_heatmapbroad %>%
  mutate(across(everything(), ~replace_na(.x, "No Rank")))
virus_heatmapbroad <- virus_heatmapbroad %>%
  mutate(
    formatted_label = paste0(subgroup, "\n", (comma(count)))
  )

## slice to top 8
virus_heatmapbroad <- virus_heatmapbroad %>% arrange(desc(count)) %>% slice_head(n = 8)

## ordering by adding a category count
virus_heatmapbroad <- virus_heatmapbroad %>% group_by(viruscategorysimple_NTorNR) %>% mutate(catcount = sum(count))
virus_heatmapbroad <- virus_heatmapbroad %>% arrange(desc(catcount),desc(count))

## adding a order id - this will be converted to a value between 0 & 1 for distinguishing between shades of non-phage using alpha
virus_heatmapbroad$subgroupid2 <- with(virus_heatmapbroad, ave(paste(subgroup), viruscategorysimple_NTorNR, FUN = function(x) match(x, unique(x))))
virus_heatmapbroad$subgroupid2 <- as.numeric(virus_heatmapbroad$subgroupid2)
virus_heatmapbroad$subgroupid <- (100 - (virus_heatmapbroad$subgroupid2)*10)/100

## also adding color scheme that is flexible with and without adapter
## make virus_heatmapbroad$viruscategorysimple_NTorNR a factor using alluvialorder2, then change the color order...
virus_heatmapbroad$viruscategorysimple_NTorNR <- factor(virus_heatmapbroad$viruscategorysimple_NTorNR, levels = alluvialorder2, ordered = TRUE)

if ("Adapter" %in% unique(virus_heatmapbroad$viruscategorysimple_NTorNR)) {
  color_values <- c("#00BFC4", "#7CAE00", "#F8766D", "#C77CFF")
} else {
  # color_values <- c("#7CAE00", "#00BFC4", "#C77CFF")
  color_values <- c("#00BFC4", "#7CAE00", "#C77CFF")
}


NTNRcontigs_treemapbroad <- ggplot2::ggplot(virus_heatmapbroad,aes(area=count,fill=viruscategorysimple_NTorNR,label=formatted_label,subgroup=viruscategorysimple_NTorNR)) + 
  treemapify::geom_treemap(layout="squarified", aes(alpha = subgroupid)) + 
  # geom_treemap(aes(alpha = subgroupid)) +
  geom_treemap_text(place = "centre", size = 18, fontface = "italic") + 
  # geom_treemap_subgroup_border() +
  # geom_treemap_subgroup_text(place = "centre", grow = TRUE, alpha = 0.7) +
  scale_fill_manual(values = color_values) +
  theme_minimal(base_size = 14, base_family = "sans") +
  guides(alpha = "none") + # this removes the alpha legend only
  labs(title="Treemap of all viral contigs found by NT + NR searches, with contig counts", fill="Virus Category") #

## saving all plots at very end

####################################################################################

## also extract all interesting virus NT blast hits - then pull as a fasta & combine with ...interesting.fasta to run the mmseqs2 easy-cluster commands:
## we first want to filter out when NR is better...  analysis_used == NTclustered
targetNT <- allchunks_diamondnr_andblastntclustered_viruses %>% dplyr::filter(analysis_used == "NTclustered")
targetNT <- targetNT %>% select(target_NTclustered) %>% unique() %>% drop_na()
targetNT <- targetNT %>% separate_wider_delim(target_NTclustered, delim = "|", names = c("remove1", "id", "remove2"), too_many = "drop")
targetNT <- targetNT %>% select(id)
write.table(targetNT, file = paste0("list_of_targetids_forclustering.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#######################################################################################
#### outputing plots and table last
ggsave(filename = paste("taxonomy_hits_viruses_alluvialplot_",Sys.Date(),".png", sep=""), alluvial_plotv, width = 18, height = 9, units = "in", limitsize = FALSE)
ggsave(filename = paste("taxonomy_hits_viruses_alluvialplot_",Sys.Date(),".pdf", sep=""), alluvial_plotv, width = 18, height = 9, units = "in", limitsize = FALSE)

write.table(data_alluvialv, file = paste0("taxonomy_hits_viruses_alluvialplot_counts_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

ggsave(filename = paste("taxonomy_hits_virus_treemap_",Sys.Date(),".png", sep=""), NTNRcontigs_treemapbroad, width = 18, height = 9, units = "in", limitsize = FALSE)
ggsave(filename = paste("taxonomy_hits_virus_treemap_",Sys.Date(),".pdf", sep=""), NTNRcontigs_treemapbroad, width = 18, height = 9, units = "in", limitsize = FALSE)

ggsave(filename = paste("taxonomy_hits_viruses_nonphage_LCA_treemap_",Sys.Date(),".png", sep=""), virus_treemap, width = 18, height = 9, units = "in", limitsize = FALSE)
ggsave(filename = paste("taxonomy_hits_viruses_nonphage_LCA_treemap_",Sys.Date(),".pdf", sep=""), virus_treemap, width = 18, height = 9, units = "in", limitsize = FALSE)

write.table(virus_heatmap, file = paste0("taxonomy_hits_viruses_nonphage_contigs_treemap_counts_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
