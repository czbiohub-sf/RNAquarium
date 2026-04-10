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
outpath <- str_c(workingpath, "/RNAquarium_outputs")
setwd(outpath)


## reload datasets from part 1
## load 'most recent version'
allchunks_diamondnr_andblastntclustered <- read_tsv("taxonomy_hits_nonhost_mostrecent.tsv.gz")


## also ALL BROAD CATEGORIES
allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "Viruses")

allchunks_diamondnr_andblastntclustered_bacteria <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "Bacteria")
allchunks_diamondnr_andblastntclustered_arthropoda <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "Arthropoda")
allchunks_diamondnr_andblastntclustered_plants <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "Plants")
allchunks_diamondnr_andblastntclustered_chordates <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "Chordata")
allchunks_diamondnr_andblastntclustered_fungi <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "Fungi")
allchunks_diamondnr_andblastntclustered_otherEukaryota <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "other_Eukaryota")
allchunks_diamondnr_andblastntclustered_SAR_Eukaryotes <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "SAR_Eukaryotes")
allchunks_diamondnr_andblastntclustered_archaea <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "Archaea")
allchunks_diamondnr_andblastntclustered_mollusca <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "Mollusca")
allchunks_diamondnr_andblastntclustered_annelida <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "Annelida")
allchunks_diamondnr_andblastntclustered_nematoda <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "Nematoda")
allchunks_diamondnr_andblastntclustered_platyhelminthes <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "Platyhelminthes")

allchunks_diamondnr_andblastntclustered_allbutchordates <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR != "Chordata")


## seqtk commands will be in a separate script - part 2
library(phylotools)

# Check if masked file exists and load appropriate fasta
if (file.exists("taxonomy_hits_nonhost_list_masked.fasta")) {
  cat("BBDuk masked version of nonhost fasta found, using that fasta\n")
  fasta_nofishnohuman <- read.fasta("taxonomy_hits_nonhost_list_masked.fasta")
} else {
  cat("No masked nonhost fasta found, using taxonomy_hits_nonhost_list.fasta\n")
  fasta_nofishnohuman <- read.fasta("taxonomy_hits_nonhost_list.fasta")
}

## rename seq.name to fullquery
fasta_nofishnohuman <- fasta_nofishnohuman %>% rename(query = seq.name)


allchunks_diamondnr_andblastntclustered_viruses_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_viruses, fasta_nofishnohuman)
write.table(allchunks_diamondnr_andblastntclustered_viruses_withsequence, file = paste0("taxonomy_hits_viruses0_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

allchunks_diamondnr_andblastntclustered_fungi_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_fungi, fasta_nofishnohuman)
write.table(allchunks_diamondnr_andblastntclustered_fungi_withsequence, file = paste0("taxonomy_hits_fungi_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

allchunks_diamondnr_andblastntclustered_archaea_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_archaea, fasta_nofishnohuman)
write.table(allchunks_diamondnr_andblastntclustered_archaea_withsequence, file = paste0("taxonomy_hits_archaea_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

allchunks_diamondnr_andblastntclustered_platyhelminthes_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_platyhelminthes, fasta_nofishnohuman)
write.table(allchunks_diamondnr_andblastntclustered_platyhelminthes_withsequence, file = paste0("taxonomy_hits_platyhelminthes_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

allchunks_diamondnr_andblastntclustered_mollusca_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_mollusca, fasta_nofishnohuman)
write.table(allchunks_diamondnr_andblastntclustered_mollusca_withsequence, file = paste0("taxonomy_hits_mollusca_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

allchunks_diamondnr_andblastntclustered_arthropoda_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_arthropoda, fasta_nofishnohuman)
write.table(allchunks_diamondnr_andblastntclustered_arthropoda_withsequence, file = paste0("taxonomy_hits_arthropoda_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

allchunks_diamondnr_andblastntclustered_nematoda_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_nematoda, fasta_nofishnohuman)
write.table(allchunks_diamondnr_andblastntclustered_nematoda_withsequence, file = paste0("taxonomy_hits_nematoda_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

allchunks_diamondnr_andblastntclustered_annelida_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_annelida, fasta_nofishnohuman)
write.table(allchunks_diamondnr_andblastntclustered_annelida_withsequence, file = paste0("taxonomy_hits_annelida_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

allchunks_diamondnr_andblastntclustered_chordates_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_chordates, fasta_nofishnohuman)
write.table(allchunks_diamondnr_andblastntclustered_chordates_withsequence, file = paste0("taxonomy_hits_chordates_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

allchunks_diamondnr_andblastntclustered_bacteria_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_bacteria, fasta_nofishnohuman)
write.table(allchunks_diamondnr_andblastntclustered_bacteria_withsequence, file = paste0("taxonomy_hits_bacteria_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

allchunks_diamondnr_andblastntclustered_plants_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_plants, fasta_nofishnohuman)
write.table(allchunks_diamondnr_andblastntclustered_plants_withsequence, file = paste0("taxonomy_hits_plants_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

allchunks_diamondnr_andblastntclustered_otherEukaryota_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_otherEukaryota, fasta_nofishnohuman)
write.table(allchunks_diamondnr_andblastntclustered_otherEukaryota_withsequence, file = paste0("taxonomy_hits_otherEukaryota_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

allchunks_diamondnr_andblastntclustered_SAR_Eukaryotes_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_SAR_Eukaryotes, fasta_nofishnohuman)
write.table(allchunks_diamondnr_andblastntclustered_SAR_Eukaryotes_withsequence, file = paste0("taxonomy_hits_SAReukaryotes_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

allchunks_diamondnr_andblastntclustered_withsequence <- left_join(allchunks_diamondnr_andblastntclustered, fasta_nofishnohuman)
write_tsv(allchunks_diamondnr_andblastntclustered_withsequence, file = paste0("taxonomy_hits_nonhost_withsequence_mostrecent.tsv.gz"))

### adding no chordates
allchunks_diamondnr_andblastntclustered_allbutchordates_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_allbutchordates, fasta_nofishnohuman)
write_tsv(allchunks_diamondnr_andblastntclustered_allbutchordates_withsequence, file = paste0("taxonomy_hits_nonhostnochordates_withsequence_mostrecent.tsv.gz"))

## checking function

# Generic function to check missing values in the rightmost column
check_join_success <- function(df, df_name = "dataframe") {
  # Get the rightmost column
  rightmost_col <- ncol(df)
  col_name <- names(df)[rightmost_col]
  
  # Count missing values
  missing_count <- df %>%
    summarise(missing = sum(is.na(.data[[col_name]]))) %>%
    pull(missing)
  
  # Print results
  cat("DataFrame:", df_name, "\n")
  cat("  Rightmost column:", col_name, "\n")
  cat("  Missing values:", missing_count, "out of", nrow(df), "rows\n")
  cat("  Success rate:", round((nrow(df) - missing_count) / nrow(df) * 100, 1), "%\n\n")
  
  return(missing_count)
}


# Automatically find and check all dataframes ending in "_withsequence"
withsequence_dfs <- ls(pattern = "_withsequence$")

if (length(withsequence_dfs) > 0) {
  cat("Checking", length(withsequence_dfs), "dataframes ending in '_withsequence':\n\n")
  
  # Function to extract last two parts of dataframe name
  extract_short_name <- function(full_name) {
    parts <- str_split(full_name, "_")[[1]]
    if (length(parts) >= 2) {
      return(paste(parts[(length(parts)-1):length(parts)], collapse = "_"))
    } else {
      return(full_name)
    }
  }
  
  # Run check_join_success on each dataframe
  for (df_name in withsequence_dfs) {
    df_obj <- get(df_name)  # Get the actual dataframe object
    short_name <- extract_short_name(df_name)
    check_join_success(df_obj, short_name)
  }
} else {
  cat("No dataframes ending in '_withsequence' found in environment.\n")
}
