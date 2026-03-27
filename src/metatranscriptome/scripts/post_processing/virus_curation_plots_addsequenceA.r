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

library(taxonomizr)
# Define path to your taxonomizr database
db_path <- "//hpc/scratch/group.data.science/eric_temp/databases/taxonomizr/aug2025taxonomy/nameNode.sqlite"

# Step 1a: Get unique taxids from taxid_NTclustered
unique_taxids_NTclustered <- allchunks_diamondnr_andblastntclustered %>%
  filter(!is.na(taxid_NTclustered)) %>%
  distinct(taxid_NTclustered) %>%
  pull(taxid_NTclustered)

# Step 1b: Get unique taxids from taxid_NR
unique_taxids_NR <- allchunks_diamondnr_andblastntclustered %>%
  filter(!is.na(taxid_NR)) %>%
  distinct(taxid_NR) %>%
  pull(taxid_NR)

cat("Unique taxids found:\n")
cat("NTclustered:", length(unique_taxids_NTclustered), "\n")
cat("NR:", length(unique_taxids_NR), "\n")

# Step 2a: Get taxonomy information for NTclustered taxids
taxonomy_lookup_NTclustered <- map_dfr(unique_taxids_NTclustered, function(taxid) {
  tax_info <- tryCatch({
    getTaxonomy(taxid, db_path, desiredTaxa = c("acellular root", "realm"))
  }, error = function(e) {
    data.frame(`acellular root` = NA, realm = NA, check.names = FALSE)
  })
  
  if (is.matrix(tax_info) || is.data.frame(tax_info)) {
    tax_df <- as.data.frame(tax_info, stringsAsFactors = FALSE)
    tax_df$taxid_NTclustered <- taxid
    return(tax_df)
  } else {
    return(data.frame(taxid_NTclustered = taxid, 
                      `acellular root` = NA, 
                      realm = NA, 
                      check.names = FALSE))
  }
}) %>%
  select(taxid_NTclustered, realm) %>%
  rename(tax_realm_NTclustered = realm)

# Step 2b: Get taxonomy information for NR taxids
taxonomy_lookup_NR <- map_dfr(unique_taxids_NR, function(taxid) {
  tax_info <- tryCatch({
    getTaxonomy(taxid, db_path, desiredTaxa = c("acellular root", "realm"))
  }, error = function(e) {
    data.frame(`acellular root` = NA, realm = NA, check.names = FALSE)
  })
  
  if (is.matrix(tax_info) || is.data.frame(tax_info)) {
    tax_df <- as.data.frame(tax_info, stringsAsFactors = FALSE)
    tax_df$taxid_NR <- taxid
    return(tax_df)
  } else {
    return(data.frame(taxid_NR = taxid, 
                      `acellular root` = NA, 
                      realm = NA, 
                      check.names = FALSE))
  }
}) %>%
  select(taxid_NR, realm) %>%
  rename(tax_realm_NR = realm)

# Step 3: Join both realm columns back to your main dataframe
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>%
  left_join(taxonomy_lookup_NTclustered, by = "taxid_NTclustered") %>%
  left_join(taxonomy_lookup_NR, by = "taxid_NR")

# Step 4a: Replace NAs in tax_clade_NTclustered with values from tax_realm_NTclustered
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>%
  mutate(tax_clade_NTclustered = ifelse(is.na(tax_clade_NTclustered), 
                                        tax_realm_NTclustered, 
                                        tax_clade_NTclustered))

# Step 4b: Replace NAs in tax_clade_NR with values from tax_realm_NR
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>%
  mutate(tax_clade_NR = ifelse(is.na(tax_clade_NR), 
                               tax_realm_NR, 
                               tax_clade_NR))


allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(tax_clade_NR, .after = tax_superkingdom_NR)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(tax_clade_NTclustered, .after = tax_superkingdom_NTclustered)

## finally now remove tax_realm_NTclustered & tax_realm_NR
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% select(-tax_realm_NTclustered) %>% select(-tax_realm_NR)

## now recreate tax_realm_NTorNR
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>%
  mutate(
    tax_clade_NTorNR = if_else(analysis_used == "NR", tax_clade_NR, tax_clade_NTclustered)
  )

## relocate after tax_superkingdom_NTorNR
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(tax_clade_NTorNR, .after = tax_superkingdom_NTorNR)

## then save only files before removing columns here, then rest at end
write.table(allchunks_diamondnr_andblastntclustered, file = paste0("taxonomy_hits_viruses0_fullcols_mostrecent2.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
