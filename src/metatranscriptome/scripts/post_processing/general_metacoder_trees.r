library(data.table)
library(scales)
library(fs)
library(taxize)
library(taxonomizr)
library(tidyverse)
###
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

library(phylotools)
library(metacoder)


workingpath <- getwd()
workingpathdash <- str_c(workingpath, "/")

## saving all outputs to RNAquarium_outputs/metacoder_trees
outpath <- str_c(workingpath, "/RNAquarium_outputs")
#setwd(outpath)
outpathmetacodertrees <- str_c(outpath, "/metacoder_trees")
outpathvirus <- str_c(outpath, "/virus_outputs")

## set wd to RNAquarium_outputs
setwd(outpath)

# Set the folder path
#folder_path <- "path/to/your/folder"

# Get all .tsv files (but not .tsv.gz)
tsv_files <- list.files(outpath, pattern = "\\.tsv$", full.names = TRUE)

# Filter out:
# - .tsv.gz files
# - taxonomy_hits_viruses0 file also any taxonomy_hits_nonhost_
## also any with sequence
tsv_files <- tsv_files[!grepl("\\.tsv\\.gz$", tsv_files)]
tsv_files <- tsv_files[!grepl("taxonomy_hits_viruses0_", basename(tsv_files))]
tsv_files <- tsv_files[!grepl("taxonomy_hits_nonhost_", basename(tsv_files))]
tsv_files <- tsv_files[!grepl("_withsequence_", basename(tsv_files))]

## ALSO GRAB taxonomy_hits_viruses_mostrecent.tsv FROM virus_outputs
#setwd(outpathvirus)
tsv_files2 <- list.files(outpathvirus, pattern = "taxonomy_hits_viruses_mostrecent.tsv", full.names = TRUE)
tsv_filesall <- c(tsv_files,tsv_files2)

cat("Found", length(tsv_filesall), "TSV files\n")
tsv_filesall

## now set output to new metacoder directory
setwd(outpathmetacodertrees)

## NOW LOOP

# Loop through each file - should be tsv_filesall, not just tsv_files!
for (file_path in tsv_filesall) {
  # Get filename without path
  filename <- basename(file_path)
  
  # Extract the third part (e.g., "fungi", "arthropoda")
  # Pattern: taxonomy_hits_NAME_DATE.tsv
  parts <- str_split(filename, "_")[[1]]
  dataset_name <- parts[3]  # Third part
  
  cat("\n========================================\n")
  cat("=== Processing:", dataset_name, "===\n")
  cat("File:", filename, "\n")
  cat("========================================\n")
  
  # Read the data
  current_data <- read_tsv(file_path, show_col_types = FALSE)
  
  ####################################
  ## DATA PREPROCESSING
  ####################################
  
  # Remove clade info for bacteria and archaea
  if (dataset_name %in% c("bacteria", "archaea", "plants")) {
    cat("  Removing tax_clade_NTorNR for", dataset_name, "\n")
    current_data$tax_clade_NTorNR <- NA
  }

  # Filter out specific phyla for otherEukaryota
  if (dataset_name == "otherEukaryota") {
    cat("  Filtering out Annelida, Arthropoda, Chordata, Mollusca, Platyhelminthes, Nematoda\n")
  
    terms_to_remove <- c("Annelida", "Arthropoda", "Chordata", "Mollusca", 
                         "Platyhelminthes", "Nematoda")
  
  # Get initial row count
    initial_rows <- nrow(current_data)
    
  # Create regex pattern
    pattern <- paste(terms_to_remove, collapse = "|")
    
  # Filter: keep rows that DON'T contain any of these terms in any tax_* column
    current_data <- current_data %>%
      filter(!if_any(
          starts_with("tax_"), 
          ~str_detect((as.character(.x)), pattern)))
  
    cat("  Removed", initial_rows - nrow(current_data), "rows\n")
    cat("  Remaining rows:", nrow(current_data), "\n")
  }
  

    # Filter out rows with "Candidatus" in any taxonomy column - moving to only in bacteria
  if (dataset_name == "bacteria") {
      cat("  Filtering out rows with 'Candidatus' in taxonomy...\n")
      rows_before_candidatus <- nrow(current_data)
  
      current_data <- current_data %>%
        filter(
          !if_any(
            c(tax_superkingdom_NTorNR, tax_clade_NTorNR, tax_kingdom_NTorNR, 
              tax_phylum_NTorNR, tax_class_NTorNR, tax_order_NTorNR, 
              tax_family_NTorNR, tax_genus_NTorNR),
            ~!is.na(.x) & str_detect(as.character(.x), "Candidatus")
          )
        )   
  
      cat("  Removed", rows_before_candidatus - nrow(current_data), "rows with 'Candidatus'\n")
  }


  if (dataset_name == "viruses") {
      cat("  Filtering out rows with whole word 'clade' or 'group' in taxonomy...\n")
      rows_before_viralfilter <- nrow(current_data)


      current_data <- current_data %>%
        filter(
          !if_any(
            c(tax_superkingdom_NTorNR, tax_clade_NTorNR, tax_kingdom_NTorNR, 
              tax_phylum_NTorNR, tax_class_NTorNR, tax_order_NTorNR, 
              tax_family_NTorNR),
            ~!is.na(.x) & str_detect(as.character(.x), regex("\\bclade\\b|\\bgroup\\b", ignore_case = TRUE))
          )
        ) %>%
#              Filter out rows containing these phrases in tax_clade_NTorNR
        filter(
          is.na(tax_clade_NTorNR) | 
          !str_detect(tax_clade_NTorNR, 
                      "Norovirus GI|Norovirus GII|Human|recombinant Vesiculovirus|papillomaviruses|Grandeviruses|amphibian retroviruses|recombinants")
        )
      
      cat("  Removed", rows_before_viralfilter - nrow(current_data), "rows with problematic taxonomy\n")
  }    
    
#   # Filter out rows with "Candidatus" in any taxonomy column - moving to only in bacteria
#   cat("  Filtering out rows with 'Candidatus' in taxonomy...\n")
#   rows_before_candidatus <- nrow(current_data)
  
#   current_data <- current_data %>%
#     filter(
#       !if_any(
#         c(tax_superkingdom_NTorNR, tax_clade_NTorNR, tax_kingdom_NTorNR, 
#           tax_phylum_NTorNR, tax_class_NTorNR, tax_order_NTorNR, 
#           tax_family_NTorNR, tax_genus_NTorNR),
#         ~!is.na(.x) & str_detect(as.character(.x), "Candidatus")
#       )
#     )
  
#   cat("  Removed", rows_before_candidatus - nrow(current_data), "rows with 'Candidatus'\n")

  # Filter out rows with NA order
  cat("  Filtering out rows with NA order...\n")
  rows_before_na <- nrow(current_data)

  current_data <- current_data %>%
    filter(!is.na(tax_order_NTorNR))

  cat("  Removed", rows_before_na - nrow(current_data), "rows with NA order\n")
  cat("  Remaining rows:", nrow(current_data), "\n")

    
  current_data <- current_data %>% arrange(desc(rel_abundance))
  current_data$taxname_lca_NTorNR <- gsub("\\(", " ", current_data$taxname_lca_NTorNR)
  current_data$taxname_lca_NTorNR <- gsub("\\)", " ", current_data$taxname_lca_NTorNR)
  
  # Define columns to remove
  cols_to_remove <- c(
    "node", "coverage", "taxid_NTorNR", "evalue_NTorNR", 
    "bits_NTclustered", "bits_NR", "evalue_NTclustered", "evalue_NR",
    "target_NTclustered", "gene_NTorNR", "allele_NTorNR", "pident_NTorNR",
    "sumperc_cov_NTorNR", "alnlen_NTorNR", "mismatch_NTorNR", "qcov_NTorNR",
    "gapopen_NTorNR", "qstart_NTorNR", "qend_NTorNR", "tstart_NTorNR",
    "tend_NTorNR", "target_title_NTorNR", "sumalnlen_NTorNR", 
    "bits_percmax_NTorNR", "target_NTorNR"
  )
  
  # Remove columns (only if they exist)
  current_data <- current_data %>%
    select(-any_of(cols_to_remove))
  
  ####################################
  ## CREATING METACODER OBJECTS
  ####################################
  
  # Subset c) Order level and up - use all data, just stop at order
  cat("\nCreating order-and-up subset...\n")
  
  # Create classification strings UP TO ORDER for ALL rows
  temp_data_order <- current_data %>%
    rowwise() %>%
    mutate(
      classification = paste(
        na.omit(c(
          tax_superkingdom_NTorNR, 
          tax_clade_NTorNR,
          tax_kingdom_NTorNR,
          tax_phylum_NTorNR,
          tax_class_NTorNR,
          tax_order_NTorNR
        )),
        collapse = ";"
      )
    ) %>%
    ungroup()
  
  cat("  Total rows:", nrow(temp_data_order), "\n")
  
  # Get unique classifications
  unique_classifications_order <- temp_data_order %>%
    distinct(classification) %>%
    filter(classification != "")
  
  cat("  Unique order-level classifications:", nrow(unique_classifications_order), "\n")
  
  # Check depth distribution
  unique_classifications_order %>%
    mutate(depth = str_count(classification, ";") + 1) %>%
    count(depth) %>%
    print()
  
  # Parse to create taxmap
  obj_order_up <- unique_classifications_order %>%
    parse_tax_data(
      class_cols = "classification",
      class_sep = ";"
    )
  
  cat("  Parsed tree has", length(obj_order_up$taxon_ids()), "taxa\n")
  
  # Add full data as query_data
  obj_order_up$data$query_data <- temp_data_order %>%
    left_join(
      obj_order_up$data$tax_data %>% select(taxon_id, classification),
      by = "classification"
    ) %>%
    filter(!is.na(taxon_id))
  
  cat("  Query data has", nrow(obj_order_up$data$query_data), "rows\n")
  
  # Verify depth distribution
  cat("\nFinal depth distribution:\n")
  obj_order_up$data$tax_data %>%
    mutate(depth = str_count(classification, ";") + 1) %>%
    count(depth) %>%
    print()
  
  cat("  Filtered to", length(obj_order_up$taxon_ids()), "taxa\n")
  
  # Verify - Check leaf names at each depth level
  leaf_check <- obj_order_up$data$tax_data %>%
    filter(!taxon_id %in% obj_order_up$edge_list$from) %>%
    mutate(
      depth = str_count(classification, ";") + 1,
      name = sapply(strsplit(classification, ";"), function(x) tail(x, 1))
    )
  
  cat("\nLeaf node depth distribution:\n")
  print(count(leaf_check, depth))
  
  cat("\nSample leaf names by depth:\n")
  for (d in 1:6) {
    if (any(leaf_check$depth == d)) {
      cat("\n--- Depth", d, "---\n")
      leaf_check %>%
        filter(depth == d) %>%
        select(name) %>%
        head(10) %>%
        print()
    }
  }
  
  ## Now calculating abundances from object
  cat("\n=== Calculating abundances ===\n")
  
  # Identify all numeric columns
  numeric_cols <- names(obj_order_up$data$query_data)[sapply(obj_order_up$data$query_data, is.numeric)]
  cat("Numeric columns found:", paste(numeric_cols, collapse = ", "), "\n")
  
  # Calculate for order and up
  cat("Calculating abundances for filtered to order subset...\n")
  obj_order_up$data$tax_abund <- calc_taxon_abund(
    obj_order_up,
    "query_data",
    cols = numeric_cols
  )
  
  cat("\n=== Setup complete! ===\n")
  cat("  - obj_order_up: Order level and up with", length(obj_order_up$taxon_ids()), "taxa\n")
  
  ####################################
  ## METACODER HEAT TREES
  ####################################
  
  # Set intervals based on dataset
  if (dataset_name %in% c("chordates", "plants")) {
    interval_range <- c(1e1, 1e8)
    cat("  Using extended interval range (1e1 to 1e8) for", dataset_name, "\n")
  } else {
    interval_range <- c(1e1, 1e7)
  }
        
  #### FIRST BY BITSCORE
  cat("\n=== Creating heat trees by bitscore ===\n")
  set.seed(11)
  
  heattrees_order <- heat_tree(
    obj_order_up,
    node_size = obj_order_up$data$tax_abund$bits_NTorNR,
    node_color = obj_order_up$data$tax_abund$bits_NTorNR,
    node_label = taxon_names,
    node_label_size_range = c(0.005, 0.015),
    node_label_max = 800,
    tree_label_max = 800,
    node_color_digits = 0,
    node_size_digits = 0,
#    edge_size_range = c(0.001, 0.03),  # <-- Add for thin branches? was 0.003 or try set just min or max leaving NA for other
    repel_force = 1.5,
    node_size_interval = interval_range,
    node_color_interval = interval_range,
    node_size_axis_label = "Bitscore",
    node_color_axis_label = "Bitscore"
  )
  
  heattrees_orderg <- heat_tree(
    obj_order_up,
    node_size = obj_order_up$data$tax_abund$bits_NTorNR,
    node_color = obj_order_up$data$tax_abund$bits_NTorNR,
    node_label = taxon_names,
    node_label_size_range = c(0.005, 0.015),
    node_label_max = 800,
    tree_label_max = 800,
    node_color_digits = 0,
    node_size_digits = 0,
#    edge_size_range = c(0.001, 0.03),  # <-- Add this for thin branches
    repel_force = 1.5,
    node_size_interval = interval_range,
    node_color_interval = interval_range,
    background_color = "gray95", 
    node_size_axis_label = "Bitscore",
    node_color_axis_label = "Bitscore"
  )
  
  heattrees_orderDH <- heat_tree(
    obj_order_up,
    node_size = obj_order_up$data$tax_abund$bits_NTorNR,
    node_color = obj_order_up$data$tax_abund$bits_NTorNR,
    node_label = taxon_names,
    node_label_size_range = c(0.005, 0.015),
    node_label_max = 800,
    tree_label_max = 800,
    node_color_digits = 0,
    node_size_digits = 0,
    repel_force = 1.5,
#    edge_size_range = c(0.001, 0.03),  # <-- Add this for thin branches
    node_size_interval = interval_range,
    node_color_interval = interval_range,
    layout = "davidson-harel",
    initial_layout = "reingold-tilford",
    node_size_axis_label = "Bitscore",
    node_color_axis_label = "Bitscore"
  )
  
  heattrees_orderDHg <- heat_tree(
    obj_order_up,
    node_size = obj_order_up$data$tax_abund$bits_NTorNR,
    node_color = obj_order_up$data$tax_abund$bits_NTorNR,
    node_label = taxon_names,
    node_label_size_range = c(0.005, 0.015),
    node_label_max = 800,
    tree_label_max = 800,
    node_color_digits = 0,
    node_size_digits = 0,
#    edge_size_range = c(0.001, 0.03),  # <-- Add this for thin branches
    repel_force = 1.5,
    node_size_interval = interval_range,
    node_color_interval = interval_range,
    background_color = "gray95", 
    layout = "davidson-harel",
    initial_layout = "reingold-tilford",
    node_size_axis_label = "Bitscore",
    node_color_axis_label = "Bitscore"
  )
  
  # Save bitscore plots
  ggsave(filename = paste0("taxonomy_hits_", dataset_name, "_bitscore_order_", Sys.Date(), ".png"), 
         heattrees_order, width = 28, height = 25, units = "in", limitsize = FALSE)
  ggsave(filename = paste0("taxonomy_hits_", dataset_name, "_bitscore_order_davidson-harel_layout_", Sys.Date(), ".png"), 
         heattrees_orderDH, width = 28, height = 25, units = "in", limitsize = FALSE)
  ggsave(filename = paste0("taxonomy_hits_", dataset_name, "_bitscore_ordergray_", Sys.Date(), ".pdf"), 
         heattrees_orderg, width = 28, height = 25, units = "in", limitsize = FALSE)
  ggsave(filename = paste0("taxonomy_hits_", dataset_name, "_bitscore_ordergray_davidson-harel_layout_", Sys.Date(), ".pdf"), 
         heattrees_orderDHg, width = 28, height = 25, units = "in", limitsize = FALSE)
  
  #### THEN BY LENGTH
  cat("\n=== Creating heat trees by length ===\n")
  set.seed(11)
  
  heattreesl_order <- heat_tree(
    obj_order_up,
    node_size = obj_order_up$data$tax_abund$length,
    node_color = obj_order_up$data$tax_abund$length,
    node_label = taxon_names,
    node_label_size_range = c(0.005, 0.015),
    node_label_max = 800,
    tree_label_max = 800,
    node_color_digits = 0,
    node_size_digits = 0,
#    edge_size_range = c(0.001, 0.03),  # <-- Add this for thin branches
    repel_force = 1.5,
    node_size_interval = interval_range,
    node_color_interval = interval_range,
    node_size_axis_label = "Total length",
    node_color_axis_label = "Total length"
  )
  
  heattreesl_orderg <- heat_tree(
    obj_order_up,
    node_size = obj_order_up$data$tax_abund$length,
    node_color = obj_order_up$data$tax_abund$length,
    node_label = taxon_names,
    node_label_size_range = c(0.005, 0.015),
    node_label_max = 800,
    tree_label_max = 800,
    node_color_digits = 0,
    node_size_digits = 0,
#    edge_size_range = c(0.001, 0.03),  # <-- Add this for thin branches
    repel_force = 1.5,
    node_size_interval = interval_range,
    node_color_interval = interval_range,
    background_color = "gray95", 
    node_size_axis_label = "Total length",
    node_color_axis_label = "Total length"
  )
  
  heattreesl_orderDH <- heat_tree(
    obj_order_up,
    node_size = obj_order_up$data$tax_abund$length,
    node_color = obj_order_up$data$tax_abund$length,
    node_label = taxon_names,
    node_label_size_range = c(0.005, 0.015),
    node_label_max = 800,
    tree_label_max = 800,
    node_color_digits = 0,
    node_size_digits = 0,
#    edge_size_range = c(0.001, 0.03),  # <-- Add this for thin branches
    repel_force = 1.5,
    node_size_interval = interval_range,
    node_color_interval = interval_range,
    layout = "davidson-harel",
    initial_layout = "reingold-tilford",
    node_size_axis_label = "Total length",
    node_color_axis_label = "Total length"
  )
  
  heattreesl_orderDHg <- heat_tree(
    obj_order_up,
    node_size = obj_order_up$data$tax_abund$length,
    node_color = obj_order_up$data$tax_abund$length,
    node_label = taxon_names,
    node_label_size_range = c(0.005, 0.015),
    node_label_max = 800,
    tree_label_max = 800,
    node_color_digits = 0,
    node_size_digits = 0,
#    edge_size_range = c(0.001, 0.03),  # <-- Add this for thin branches
    repel_force = 1.5,
    node_size_interval = interval_range,
    node_color_interval = interval_range,
    background_color = "gray95", 
    layout = "davidson-harel",
    initial_layout = "reingold-tilford",
    node_size_axis_label = "Total length",
    node_color_axis_label = "Total length"
  )
  
  # Save length plots
  ggsave(filename = paste0("taxonomy_hits_", dataset_name, "_length_order_", Sys.Date(), ".png"), 
         heattreesl_order, width = 28, height = 25, units = "in", limitsize = FALSE)
  ggsave(filename = paste0("taxonomy_hits_", dataset_name, "_length_order_davidson-harel_layout_", Sys.Date(), ".png"), 
         heattreesl_orderDH, width = 28, height = 25, units = "in", limitsize = FALSE)
  ggsave(filename = paste0("taxonomy_hits_", dataset_name, "_length_ordergray_", Sys.Date(), ".pdf"), 
         heattreesl_orderg, width = 28, height = 25, units = "in", limitsize = FALSE)
  ggsave(filename = paste0("taxonomy_hits_", dataset_name, "_length_ordergray_davidson-harel_layout_", Sys.Date(), ".pdf"), 
         heattreesl_orderDHg, width = 28, height = 25, units = "in", limitsize = FALSE)
  
  #### THEN BY ABUNDANCE
  cat("\n=== Creating heat trees by abundance ===\n")
  set.seed(11)
  
  heattreea_order <- heat_tree(
    obj_order_up,
    node_size = obj_order_up$data$tax_abund$rel_abundance,
    node_color = obj_order_up$data$tax_abund$rel_abundance,
    node_label = taxon_names,
    node_label_size_range = c(0.005, 0.015),
    node_label_max = 800,
    tree_label_max = 800,
    node_color_digits = 0,
    node_size_digits = 0,
#    edge_size_range = c(0.001, 0.03),  # <-- Add this for thin branches
    repel_force = 1.5,
    node_size_interval = interval_range,
    node_color_interval = interval_range,
    node_size_axis_label = "Relative abundance",
    node_color_axis_label = "Relative abundance"
  )
  
  heattreea_orderg <- heat_tree(
    obj_order_up,
    node_size = obj_order_up$data$tax_abund$rel_abundance,
    node_color = obj_order_up$data$tax_abund$rel_abundance,
    node_label = taxon_names,
    node_label_size_range = c(0.005, 0.015),
    node_label_max = 800,
    tree_label_max = 800,
    node_color_digits = 0,
    node_size_digits = 0,
#    edge_size_range = c(0.001, 0.03),  # <-- Add this for thin branches
    repel_force = 1.5,
    node_size_interval = interval_range,
    node_color_interval = interval_range,
    background_color = "gray95", 
    node_size_axis_label = "Relative abundance",
    node_color_axis_label = "Relative abundance"
  )
  
  heattreea_orderDH <- heat_tree(
    obj_order_up,
    node_size = obj_order_up$data$tax_abund$rel_abundance,
    node_color = obj_order_up$data$tax_abund$rel_abundance,
    node_label = taxon_names,
    node_label_size_range = c(0.005, 0.015),
    node_label_max = 800,
    tree_label_max = 800,
    node_color_digits = 0,
    node_size_digits = 0,
#    edge_size_range = c(0.001, 0.03),  # <-- Add this for thin branches
    repel_force = 1.5,
    node_size_interval = interval_range,
    node_color_interval = interval_range,
    layout = "davidson-harel",
    initial_layout = "reingold-tilford",
    node_size_axis_label = "Relative abundance",
    node_color_axis_label = "Relative abundance"
  )
  
  heattreea_orderDHg <- heat_tree(
    obj_order_up,
    node_size = obj_order_up$data$tax_abund$rel_abundance,
    node_color = obj_order_up$data$tax_abund$rel_abundance,
    node_label = taxon_names,
    node_label_size_range = c(0.005, 0.015),
    node_label_max = 800,
    tree_label_max = 800,
    node_color_digits = 0,
    node_size_digits = 0,
#    edge_size_range = c(0.001, 0.03),  # <-- Add this for thin branches
    repel_force = 1.5,
    node_size_interval = interval_range,
    node_color_interval = interval_range,
    background_color = "gray95", 
    layout = "davidson-harel",
    initial_layout = "reingold-tilford",
    node_size_axis_label = "Relative abundance",
    node_color_axis_label = "Relative abundance"
  )
  
  # Save abundance plots
  ggsave(filename = paste0("taxonomy_hits_", dataset_name, "_abundance_order_", Sys.Date(), ".png"), 
         heattreea_order, width = 28, height = 25, units = "in", limitsize = FALSE)
  ggsave(filename = paste0("taxonomy_hits_", dataset_name, "_abundance_order_davidson-harel_layout_", Sys.Date(), ".png"), 
         heattreea_orderDH, width = 28, height = 25, units = "in", limitsize = FALSE)
  ggsave(filename = paste0("taxonomy_hits_", dataset_name, "_abundance_ordergray_", Sys.Date(), ".pdf"), 
         heattreea_orderg, width = 28, height = 25, units = "in", limitsize = FALSE)
  ggsave(filename = paste0("taxonomy_hits_", dataset_name, "_abundance_ordergray_davidson-harel_layout_", Sys.Date(), ".pdf"), 
         heattreea_orderDHg, width = 28, height = 25, units = "in", limitsize = FALSE)
  
  cat("\n=== Finished processing", dataset_name, "===\n")
  cat("Created 12 heat tree plots\n")
}

cat("\n========================================\n")
cat("All files processed successfully!\n")
cat("========================================\n")
        
