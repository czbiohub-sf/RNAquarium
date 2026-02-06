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

library(phylotools)
library(metacoder)
options(ENTREZ_KEY='94ef22f4482d7bc124b52e1b9036688eb608')


workingpath <- getwd()
workingpathdash <- str_c(workingpath, "/")
diamondpath <- str_c(workingpath, "/nr_diamond")
blastpath <- str_c(workingpath, "/nt_blast")

## saving all outputs to RNAquarium_outputs/virus_outputs
outpath <- str_c(workingpath, "/RNAquarium_outputs")
#setwd(outpath)

outpathvirus <- str_c(outpath, "/virus_outputs")
setwd(outpathvirus)

####################################
## THEN VIRUS SPECIFIC CATEGORIES NEXT


## load 'most recent version'
allchunks_diamondnr_andblastntclustered_interesting <- read_tsv("taxonomy_hits_viruses_nonphage_mostrecent.tsv")


## subset to remove phage, then sort by evalue to limit Metacoder tree to 400 unique taxon names...
# comparisons_viruses <- comparisons_viruses %>% dplyr::filter(!grepl('phage|Phage', taxname_lca))
# comparisons_viruses <- comparisons_viruses %>% dplyr::filter(!grepl('Caudoviricetes', tax_class))

allchunks_diamondnr_andblastntclustered_interesting <- allchunks_diamondnr_andblastntclustered_interesting %>% dplyr::filter(viruscategorysimple_NTorNR != "Phage")


## note we are already pulling max bitscore, so this is okay to use instead of e-value cutoff!
#allchunks_diamondnr_andblastntclustered_interesting <- allchunks_diamondnr_andblastntclustered_interesting %>% arrange(desc(evalue_NTorNR))



thresholded_hit_virusesabundanttop400 <- allchunks_diamondnr_andblastntclustered_interesting %>% group_by(taxname_lca_NTorNR) %>% summarize(maxrelabundance = max(bits_NTorNR))

thresholded_hit_virusesabundanttop400 <- thresholded_hit_virusesabundanttop400 %>% arrange(desc(maxrelabundance)) %>% slice_head(n = 400)


thresholded_hit_virusesabundantfewer <- inner_join(allchunks_diamondnr_andblastntclustered_interesting, thresholded_hit_virusesabundanttop400)
n_distinct(thresholded_hit_virusesabundantfewer$taxname_lca_NTorNR)

## converting all parentheses to spaces
thresholded_hit_virusesabundantfewer$taxname_lca_NTorNR <- gsub("\\(", " ", thresholded_hit_virusesabundantfewer$taxname_lca_NTorNR)
thresholded_hit_virusesabundantfewer$taxname_lca_NTorNR <- gsub("\\)", " ", thresholded_hit_virusesabundantfewer$taxname_lca_NTorNR)

## removing more
thresholded_hit_virusesabundantfewcols <- thresholded_hit_virusesabundantfewer %>% select(-node) %>% select(-coverage) %>% select(-taxid_NTorNR) %>% select(-evalue_NTorNR) %>% select(-bits_NTclustered) %>% select(-bits_NR) %>% select(-evalue_NTclustered) %>% select(-evalue_NR) %>% select(-target_NTclustered) %>% select(-gene_NTorNR) %>% select(-allele_NTorNR) %>% select(-pident_NTorNR) %>% select(-sumperc_cov_NTorNR) %>% select(-alnlen_NTorNR) %>% select(-mismatch_NTorNR) %>% select(-qcov_NTorNR) %>% select(-gapopen_NTorNR) %>% select(-qstart_NTorNR) %>% select(-qend_NTorNR) %>% select(-tstart_NTorNR) %>% select(-tend_NTorNR) %>% select(-target_title_NTorNR) %>% select(-sumalnlen_NTorNR) %>% select(-bits_percmax_NTorNR) %>% select(-target_NTorNR)


### smaller sets in case first fails
thresholded_hit_virusesabundanttop300 <- thresholded_hit_virusesabundanttop400 %>% arrange(desc(maxrelabundance)) %>% slice_head(n = 300)
thresholded_hit_virusesabundanttop200 <- thresholded_hit_virusesabundanttop400 %>% arrange(desc(maxrelabundance)) %>% slice_head(n = 200)

thresholded_hit_virusesabundantfewer300 <- inner_join(allchunks_diamondnr_andblastntclustered_interesting, thresholded_hit_virusesabundanttop300)
thresholded_hit_virusesabundantfewer200 <- inner_join(allchunks_diamondnr_andblastntclustered_interesting, thresholded_hit_virusesabundanttop200)

thresholded_hit_virusesabundantfewer300$taxname_lca_NTorNR <- gsub("\\(", " ", thresholded_hit_virusesabundantfewer300$taxname_lca_NTorNR)
thresholded_hit_virusesabundantfewer300$taxname_lca_NTorNR <- gsub("\\)", " ", thresholded_hit_virusesabundantfewer300$taxname_lca_NTorNR)
thresholded_hit_virusesabundantfewcols300 <- thresholded_hit_virusesabundantfewer300 %>% select(-node) %>% select(-coverage) %>% select(-taxid_NTorNR) %>% select(-evalue_NTorNR) %>% select(-bits_NTclustered) %>% select(-bits_NR) %>% select(-evalue_NTclustered) %>% select(-evalue_NR) %>% select(-target_NTclustered) %>% select(-gene_NTorNR) %>% select(-allele_NTorNR) %>% select(-pident_NTorNR) %>% select(-sumperc_cov_NTorNR) %>% select(-alnlen_NTorNR) %>% select(-mismatch_NTorNR) %>% select(-qcov_NTorNR) %>% select(-gapopen_NTorNR) %>% select(-qstart_NTorNR) %>% select(-qend_NTorNR) %>% select(-tstart_NTorNR) %>% select(-tend_NTorNR) %>% select(-target_title_NTorNR) %>% select(-sumalnlen_NTorNR) %>% select(-bits_percmax_NTorNR) %>% select(-target_NTorNR)

thresholded_hit_virusesabundantfewer200$taxname_lca_NTorNR <- gsub("\\(", " ", thresholded_hit_virusesabundantfewer200$taxname_lca_NTorNR)
thresholded_hit_virusesabundantfewer200$taxname_lca_NTorNR <- gsub("\\)", " ", thresholded_hit_virusesabundantfewer200$taxname_lca_NTorNR)
thresholded_hit_virusesabundantfewcols200 <- thresholded_hit_virusesabundantfewer200 %>% select(-node) %>% select(-coverage) %>% select(-taxid_NTorNR) %>% select(-evalue_NTorNR) %>% select(-bits_NTclustered) %>% select(-bits_NR) %>% select(-evalue_NTclustered) %>% select(-evalue_NR) %>% select(-target_NTclustered) %>% select(-gene_NTorNR) %>% select(-allele_NTorNR) %>% select(-pident_NTorNR) %>% select(-sumperc_cov_NTorNR) %>% select(-alnlen_NTorNR) %>% select(-mismatch_NTorNR) %>% select(-qcov_NTorNR) %>% select(-gapopen_NTorNR) %>% select(-qstart_NTorNR) %>% select(-qend_NTorNR) %>% select(-tstart_NTorNR) %>% select(-tend_NTorNR) %>% select(-target_title_NTorNR) %>% select(-sumalnlen_NTorNR) %>% select(-bits_percmax_NTorNR) %>% select(-target_NTorNR)



n_distinct(thresholded_hit_virusesabundantfewcols$taxname_lca_NTorNR)
unique(thresholded_hit_virusesabundantfewcols$taxname_lca_NTorNR)


## now actual Metacoder code...

#obj_viraltopabundance3 <- lookup_tax_data(thresholded_hit_virusesabundantfewcols, type = "taxon_name", column = "taxname_lca_NTorNR", ask = FALSE)
## getting Error: Bad Request (HTTP 400)
# Error in curl::curl_fetch_memory(x$url$url, handle = x$url$handle) : 
#   Timeout was reached: [eutils.ncbi.nlm.nih.gov] Connection timed out after 10001 milliseconds

## put this into a try loop per https://forum.posit.co/t/trycatch-with-loop/86749/3
library(taxize)

# # Set up a retry mechanism
# max_retries <- 5  # Maximum number of retries
# retry_count <- 0  # Initialize retry counter
# success <- FALSE  # Flag to indicate success
# 
# while (!success && retry_count < max_retries) {
#   tryCatch(
#     {
#       retry_count <- retry_count + 1  # Increment retry count before the attempt
#       message(paste("Attempt", retry_count, "of", max_retries, "..."))
#       
#       # Attempt to run the command
#       obj_viraltopabundance3 <- lookup_tax_data(
#         thresholded_hit_virusesabundantfewcols,
#         type = "taxon_name",
#         column = "taxname_lca_NTorNR",
#         ask = FALSE
#       )
#       
#       # If successful
#       success <- TRUE
#       message("lookup_tax_data executed successfully.")
#     },
#     error = function(e) {
#       # Handle errors
#       message(paste("Attempt", retry_count, "failed."))
#       
#       # Check if retries remain
#       if (retry_count < max_retries) {
#         message("Retrying in 5 seconds...")
#         Sys.sleep(5)  # Wait before retrying
#       } else {
#         message("Maximum retries reached. Exiting...")
#         stop("An unexpected error occurred: ", conditionMessage(e))
#       }
#     }
#   )
# }

##########################################
### 3 sets of tries
# Define your inputs
inputs <- list(
  "thresholded_hit_virusesabundantfewcols",
  "thresholded_hit_virusesabundantfewcols300",
  "thresholded_hit_virusesabundantfewcols200"
)

max_retries <- 3  # Maximum retries per input
total_attempts <- 0  # Track total attempts
success <- FALSE  # Flag to indicate success

for (input in inputs) {
  retry_count <- 0  # Reset retry counter for each input
  
  while (!success && retry_count < max_retries) {
    total_attempts <- total_attempts + 1  # Increment total attempts
    retry_count <- retry_count + 1  # Increment retry count for the current input
    
    tryCatch(
      {
        message(paste("Attempt", total_attempts, "with input:", input, "(Retry", retry_count, "of", max_retries, ")..."))
        
        # Attempt the command with the current input
        obj_viraltopabundance3o <- lookup_tax_data(
          get(input),  # Use the current input
          type = "taxon_name",
          column = "taxname_lca_NTorNR",
          ask = FALSE
        )
        
        # If successful, set success to TRUE and break out of the loop
        success <- TRUE
        message("lookup_tax_data executed successfully.")
      },
      error = function(e) {
        # Handle errors
        message(paste("Attempt", total_attempts, "failed with input:", input, "."))
        
        if (retry_count < max_retries) {
          message("Retrying in 8 seconds...")
          Sys.sleep(8)  # Wait before retrying
        } else {
          message(paste("All", max_retries, "retries failed for input:", input, ". Moving to the next input..."))
        }
      }
    )
  }
  
  # If successful, break out of the for-loop
  if (success) break
}

# If all inputs failed
if (!success) {
  stop("Failed after 15 total attempts across all inputs.")
}

## remove "unknown taxon" if it exists?
obj_viraltopabundance3o <- obj_viraltopabundance3o %>%
  filter_ambiguous_taxa()

#obj_viraltopabundance3 <- obj_viraltopabundance3o
obj_viraltopabundance3 <- obj_viraltopabundance3o$clone(deep = TRUE)
save(obj_viraltopabundance3o, file = "obj_taxonomy_hits_viruses_nophage.Rdata")



obj_viraltopabundance3$data <- calc_taxon_abund(obj_viraltopabundance3, "query_data")
save(obj_viraltopabundance3, file = "obj_taxonomy_hits_viruses_withlengthandbitscore.Rdata")

#metatree <- heat_tree(obj_virusestest21, node_label = taxon_names, node_size = n_obs, node_color = n_obs)

## instead of abundance, both by bitscore & by length
set.seed(11) # This makes the plot appear the same each time it is run 
metatree32 <- heat_tree(obj_viraltopabundance3, 
                        node_label = taxon_names,
                        node_size = rel_abundance,  ## works after line above, note change next to rel_abundance
                        node_color = rel_abundance, 
                        node_size_axis_label = "Relative abundance",
                        node_label_size_range = c(0.005, 0.015),
                        node_label_max = 800,
                        edge_label_max = 500,
                        tree_label_max = 800,
                        node_color_digits = 0,
                        node_size_digits = 0,
                        edge_color_digits = 0,
                        edge_size_digits = 0,
                        repel_force = 1.5,
                        node_size_interval = c(1e1,1e7),
                        node_color_interval = c(1e1,1e7),
                        background_color = "gray95",
                        # node_color_axis_label = "Samples with reads",
                        layout = "davidson-harel", # The primary layout algorithm
                        initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations
ggsave(filename = paste("taxonomy_hits_viruses_nophage_metacodertree_relabundance_",Sys.Date(),".png", sep=""), metatree32, width = 44, height = 34, units = "in", limitsize = FALSE)
ggsave(filename = paste("taxonomy_hits_viruses_nophage_metacodertree_relabundance_",Sys.Date(),".pdf", sep=""), metatree32, width = 44, height = 34, units = "in", limitsize = FALSE)


metatree322 <- heat_tree(obj_viraltopabundance3, 
                         node_label = taxon_names,
                         node_size = bits_NTorNR,  ## works after line above, note change from rel_abundance
                         node_color = bits_NTorNR, 
                         node_size_axis_label = "Bitscore",
                         node_label_size_range = c(0.005, 0.015),
                         node_label_max = 800,
                         edge_label_max = 500,
                         tree_label_max = 800,
                         node_color_digits = 0,
                         node_size_digits = 0,
                         edge_color_digits = 0,
                         edge_size_digits = 0,
                         repel_force = 1.5,
                         node_size_interval = c(1e1,1e7),
                         node_color_interval = c(1e1,1e7),
                         background_color = "gray95",                         
                         # node_color_axis_label = "Samples with reads",
                         layout = "davidson-harel", # The primary layout algorithm
                         initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations
ggsave(filename = paste("taxonomy_hits_viruses_nophage_metacodertree_bitscore_",Sys.Date(),".png", sep=""), metatree322, width = 44, height = 34, units = "in", limitsize = FALSE)
ggsave(filename = paste("taxonomy_hits_viruses_nophage_metacodertree_bitscore_",Sys.Date(),".pdf", sep=""), metatree322, width = 44, height = 34, units = "in", limitsize = FALSE)

metatree323 <- heat_tree(obj_viraltopabundance3, 
                         node_label = taxon_names,
                         node_size = length,  ## works after line above, note change from rel_abundance
                         node_color = length, 
                         node_size_axis_label = "Contig length",
                         node_label_size_range = c(0.005, 0.015),
                         node_label_max = 800,
                         edge_label_max = 500,
                         tree_label_max = 800,
                         node_color_digits = 0,
                         node_size_digits = 0,
                         edge_color_digits = 0,
                         edge_size_digits = 0,
                         repel_force = 1.5,
                         node_size_interval = c(1e1,1e7),
                         node_color_interval = c(1e1,1e7),
                         background_color = "gray95",
                         # node_color_axis_label = "Samples with reads",
                         layout = "davidson-harel", # The primary layout algorithm
                         initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations
ggsave(filename = paste("taxonomy_hits_viruses_nophage_metacodertree_contiglength_",Sys.Date(),".png", sep=""), metatree323, width = 44, height = 34, units = "in", limitsize = FALSE)
ggsave(filename = paste("taxonomy_hits_viruses_nophage_metacodertree_contiglength_",Sys.Date(),".pdf", sep=""), metatree323, width = 44, height = 34, units = "in", limitsize = FALSE)

