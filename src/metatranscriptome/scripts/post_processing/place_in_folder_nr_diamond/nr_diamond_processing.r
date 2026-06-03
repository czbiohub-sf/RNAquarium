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

####UPDATE THIS WITH YOUR TAXNOMIZR PATH####  defining paths with Taxonomizr sqlite database(s)
taxonomizr_base_dir <- "/path/to/databases/taxonomizr/aug2025taxonomy"
taxonomizr_fallback_dir <- "/path/to/databases/taxonomizr/feb2025taxonomy"

taxonomizr_db_path <- file.path(taxonomizr_base_dir, "nameNode.sqlite")
taxonomizr_db_path_fallback <- file.path(taxonomizr_fallback_dir, "nameNode.sqlite")


## CHANGING RECURSE = FALSE TO PREVENT DOUBLE COUNTING FROM SUBFOLDERS
alldiamondfiles <- fs::dir_ls(workingpath, glob="*.diamond.txt.gz", recurse = FALSE)

## chunk names
chunknames <- alldiamondfiles
chunknames <- str_replace(chunknames, workingpathdash, "")
chunknames <- str_replace(chunknames, ".diamond.txt.gz", "")
chunknames <- str_replace(chunknames, "_nonzfhum.", "")

### the looped version - for a directory with many files, this will read one at a time

for(i in seq(alldiamondfiles)){
  single_hit <- fread(alldiamondfiles[i], verbose = TRUE, quote="")
  single_hit <- as_tibble(single_hit)
  
  single_hit <- single_hit %>% rename("query" = "V1") %>% 
    rename("target" = "V2") %>% 
    rename("taxid" = "V3") %>% 
    rename("taxname" = "V4") %>% 
    rename("skingdom" = "V5") %>% 
    rename("pident" = "V6") %>% 
    rename("alnlen" = "V7") %>% 
    rename("mismatch" = "V8") %>% 
    rename("qcov" = "V9") %>% 
    rename("gapopen" = "V10") %>% 
    rename("qstart" = "V11") %>% 
    rename("qend" = "V12") %>% 
    rename("tstart" = "V13") %>% 
    rename("tend" = "V14") %>% 
    rename("evalue" = "V15") %>% 
    rename("bits" = "V16") %>% 
    rename("target_title" = "V17")
  
  single_hit$analysis <- "diamond"
  
  # for Diamond NR only - check first 1000 rows for NA or "N/A" values in the taxid column
  if ("taxid" %in% names(single_hit)) {
    
    taxid_check <- single_hit %>%
      slice_head(n = 1000) %>%
      pull(taxid)
    
    has_missing <- any(is.na(taxid_check) | taxid_check == "N/A")
    
    if (has_missing) {
      stop("WARNING: Your first Diamond dataframe includes some NA values in the taxid column. Please run script slurm_diamond_taxid_update_blastdbcmd.sh in the nr_diamond folder to find taxid values from a local database. Then re-run this script.")
    } else {
      cat("TaxID check passed: no NA or 'N/A' values found in the first 1000 rows.\n")
    }
    
  } else {
    stop("ERROR: 'taxid' column not found in the 'single_hit' dataframe.")
  }
  
  # separate first to get length...
  single_hit <- single_hit %>% separate(query, into = c("bioproject", "remove", "remove1", "node", "remove2", "length", "remove3", "coverage", "gene", "allele"), sep = "_", remove = FALSE, convert = TRUE, extra = "drop")
  single_hit$remove <- NULL
  single_hit$remove1 <- NULL
  single_hit$remove2 <- NULL
  single_hit$remove3 <- NULL
  
  
  ### NOTE for nr divide alnlen by 3
  single_hit <- single_hit %>% arrange(desc(query),desc(bits))
  single_hitt <- single_hit %>% group_by(query) %>% summarize(sumalnlen = sum(alnlen * 3))
  
  single_hit <- left_join(single_hit,single_hitt)
  ## then create a new column sumperc_cov
  ## first force as.numeric
  single_hit$length <- as.numeric(single_hit$length)
  single_hit$coverage <- as.numeric(single_hit$coverage)
  
  single_hit$sumperc_cov <- single_hit$sumalnlen / single_hit$length
  single_hit <- single_hit %>% relocate(sumperc_cov, .after = pident)
  #####
  
  ## 99% maxbits filter
  single_hit_maxbits <- single_hit %>% group_by(query) %>% summarize(maxbits = max(bits))
  single_hitss <- left_join(single_hit, single_hit_maxbits)
  single_hitss$bits_percmax <- single_hitss$bits / single_hitss$maxbits
  top_hits <- single_hitss %>% dplyr::filter(bits_percmax > 0.989999)
  
  top_hits$rel_abundance <- (top_hits$length * top_hits$coverage * 0.01)
  top_hits <- top_hits %>% relocate(rel_abundance, .after = coverage)
  

  top_hits$taxid <- gsub(";.*", "", top_hits$taxid)
  top_hits$taxname <- gsub("0.*", "", top_hits$taxname)
  
  ##need to convert taxid to as.numeric
  top_hits$taxid <- as.numeric(top_hits$taxid)
  
  ## often way too many taxids, so extract only unique
  n_distinct(top_hits$taxid)
  taxid_unique <- unique(top_hits$taxid)
  taxa <- getTaxonomy(
    taxid_unique,
    taxonomizr_db_path,
    desiredTaxa = c("domain", "acellular root", "clade", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  )  
  
  taxa.df <- as.data.frame(taxa)
  taxa.df <- as_tibble(taxa, rownames = "taxid")
  taxa.df$taxid <- as.numeric(taxa.df$taxid)
  
  ## AS PART OF UPDATED TAXONOMY - NEED TO COALESCE TWO FIELDS TO RECREATE OLD SUPERKINGDOM
  #should always have one of acellular root" or "domain" as NA so thus can re-create "superkingdom" column using one simple coalesce
  ## also rename
  taxa.df <- taxa.df %>% rename(AcellularRoot = `acellular root`)
  taxa.df <- taxa.df %>% mutate(superkingdom = coalesce(domain, AcellularRoot))
  taxa.df <- taxa.df %>% select(-any_of(c("domain", "AcellularRoot")))
  ## add check - for rows with all NAs, run again using feb2025 database
  
  
  # NEW: Second pass for missing taxids
  #Find rows where all columns except taxid are NA
  all_na_rows <- taxa.df %>%
    rowwise() %>%
    mutate(all_cols_na = all(is.na(c_across(-taxid)))) %>%
    ungroup() %>%
    filter(all_cols_na == TRUE)
  
  # Extract taxids from these rows
  missing_taxids <- all_na_rows$taxid
  
  cat("Found", length(missing_taxids), "taxids with all NA values for second Taxonomizr pass\n")
  
  if (length(missing_taxids) > 0) {
    #  Second pass with older database
    taxa_second_pass <- getTaxonomy(
      missing_taxids,
      taxonomizr_db_path_fallback,
      desiredTaxa = c("superkingdom", "clade", "kingdom", "phylum", "class", "order", "family", "genus", "species")
    )
    
    taxa_second_pass.df <- as.data.frame(taxa_second_pass)
    taxa_second_pass.df <- as_tibble(taxa_second_pass, rownames = "taxid")
    taxa_second_pass.df$taxid <- as.numeric(taxa_second_pass.df$taxid)
    
    # Update the original taxa.df with results from second pass
    # Remove the all-NA rows from original dataframe
    taxa.df <- taxa.df %>%
      filter(!taxid %in% missing_taxids)
    
    # Bind the second pass results
    taxa.df <- bind_rows(taxa.df, taxa_second_pass.df)
    
    cat("Second Taxonomizr pass completed. Found taxonomy for", 
        sum(!is.na(taxa_second_pass.df$superkingdom)), 
        "out of", length(missing_taxids), "missing taxids\n")
  } else {
    cat("No taxids needing second Taxonomizr pass - all were resolved in first Taxonomizr pass\n")
  }
  
  ## now need to unlabel rows and name column
  top_hitswithtaxa <- left_join(top_hits,taxa.df)
  
  rm(top_hits)
  top_hits <- top_hitswithtaxa
  
  
  ########### END TAXONOMIZR CODE ###########
  
  ## update to now remove coalesced domain, `acellular root`
  top_hits <- top_hits %>% select(-any_of(c("skingdom", "tax_kingdom", "domain", "AcellularRoot")))
  top_hits <- top_hits %>% rename("tax_superkingdom" = "superkingdom") %>% rename("tax_clade" = "clade") %>% rename("tax_kingdom" = "kingdom") %>% rename("tax_phylum" = "phylum") %>%rename("tax_class" = "class") %>% rename("tax_order" = "order") %>% rename("tax_family" = "family") %>% rename("tax_genus" = "genus") %>% rename("tax_species" = "species")
  
  ## TWEAKED TAXONOMIZR CODE SO WE CAN SEARCH FOR SAR 
  top_hits$taxoncategory <- ifelse((grepl("Viruses", top_hits$tax_superkingdom) == TRUE), "Viruses",
                                   ifelse((grepl("Bacteria", top_hits$tax_superkingdom) == TRUE), "Bacteria", 
                                          ifelse((grepl("Fungi", top_hits$tax_kingdom) == TRUE), "Fungi", 
                                                 #ifelse((grepl("Ascomycota", top_hits$phylum) == TRUE), "Fungi",        
                                                 ifelse((grepl("Viridiplantae", top_hits$tax_kingdom) == TRUE), "Plants", 
                                                        ifelse((grepl("Sar", top_hits$tax_clade) == TRUE), "SAR_Eukaryotes",
                                                               ifelse((grepl("Actinopteri", top_hits$tax_class) == TRUE), "Fishes",
                                                                      #ifelse((grepl("Basidiomycota", top_hits$phylum) == TRUE), "Fungi",   
                                                                      ifelse((grepl("Primates", top_hits$tax_order) == TRUE), "Human", top_hits$tax_phylum)))))))
  
  
  ## then replace_na
  top_hits <- top_hits %>% mutate(across(taxoncategory, ~ replace_na(., "Root_unresolved")))
  
  ## weird one-ff if tax_kingdom is all zeros
  top_hits$tax_kingdom <- as.character(top_hits$tax_kingdom)
  top_hits <- top_hits %>%
    mutate(tax_kingdom = na_if(tax_kingdom, "0"))
  top_hits <- top_hits %>%
    mutate(tax_phylum = na_if(tax_phylum, "0"))
  
  thresholded_hit_nofishnomammals0 <- top_hits %>% filter(!grepl("Actinopteri", top_hits$tax_class))
  thresholded_hit_nofishnomammals0 <- thresholded_hit_nofishnomammals0 %>% filter(!grepl("Primates", thresholded_hit_nofishnomammals0$tax_order))
  ## ADDING VERY LAST - ACCOUNTING FOR FISH OR PRIMATE HITS FOR DARK MATTER    
  currentchunkname <- chunknames[i]
  thresholded_hit_fishorprimates1 <- top_hits %>% filter(grepl("Actinopteri", top_hits$tax_class))
  thresholded_hit_fishorprimates2 <- top_hits %>% filter(grepl("Primates", top_hits$tax_order))
  thresholded_hit_fishorprimates <- bind_rows(thresholded_hit_fishorprimates1,thresholded_hit_fishorprimates2)
  thresholded_hit_fishorprimates <- thresholded_hit_fishorprimates %>% select(query)
  thresholded_hit_fishorprimates <- thresholded_hit_fishorprimates %>% unique()
  ## then save as a txt file
  write.table(thresholded_hit_fishorprimates, file = paste0(currentchunkname,"_diamond_hits_fishorprimates",".txt"), sep = "\t", row.names = FALSE, quote = FALSE, , col.names = FALSE)
  
  #################### CODE TO PULL LCA FROM TAXONOMY ##############################
  
  ## first arrange thresholded_hit_nofishnomammals0
  thresholded_hit_nofishnomammals0 <- thresholded_hit_nofishnomammals0 %>% arrange(bioproject,desc(length),desc(node),desc(bits))
  ## now slice by query and taxid
  thresholded_hit_nofishnomammals7 <- thresholded_hit_nofishnomammals0 %>% group_by(query,taxid) %>% slice_head(n = 1) %>% ungroup()

  thresholded_hit_nofish_dups <- thresholded_hit_nofishnomammals7 %>% group_by(query) %>% add_tally() %>% dplyr::filter(n > 1) %>% ungroup()
  thresholded_hit_nofish_rest <- thresholded_hit_nofishnomammals7 %>% group_by(query) %>% add_tally() %>% dplyr::filter(n == 1) %>% ungroup()

  thresholded_hit_nofish_dups <- thresholded_hit_nofish_dups %>% dplyr::group_by(query, tax_superkingdom) %>%
    tidyr::fill(tax_superkingdom, .direction = "downup") %>%
    tidyr::fill(tax_clade, .direction = "downup") %>%
    tidyr::fill(tax_kingdom, .direction = "downup") %>%
    tidyr::fill(tax_phylum, .direction = "downup") %>%
    tidyr::fill(tax_class, .direction = "downup") %>%
    tidyr::fill(tax_order, .direction = "downup") %>%
    tidyr::fill(tax_family, .direction = "downup") %>%
    tidyr::fill(tax_genus, .direction = "downup") %>%
    dplyr::ungroup()
  
  examining_dups1 <- thresholded_hit_nofish_dups %>% group_by(query) %>% summarise(across(tax_superkingdom, n_distinct))
  examining_dups1 <- examining_dups1 %>% rename(taxlineage1_n = tax_superkingdom)
  examining_dups2 <- thresholded_hit_nofish_dups %>% group_by(query) %>% summarise(across(tax_clade, n_distinct))
  examining_dups2 <- examining_dups2 %>% rename(taxlineage2_n = tax_clade)
  examining_dups3 <- thresholded_hit_nofish_dups %>% group_by(query) %>% summarise(across(tax_kingdom, n_distinct))
  examining_dups3 <- examining_dups3 %>% rename(taxlineage3_n = tax_kingdom)
  examining_dups4 <- thresholded_hit_nofish_dups %>% group_by(query) %>% summarise(across(tax_phylum, n_distinct))
  examining_dups4 <- examining_dups4 %>% rename(taxlineage4_n = tax_phylum)
  examining_dups5 <- thresholded_hit_nofish_dups %>% group_by(query) %>% summarise(across(tax_class, n_distinct))
  examining_dups5 <- examining_dups5 %>% rename(taxlineage5_n = tax_class)
  examining_dups6 <- thresholded_hit_nofish_dups %>% group_by(query) %>% summarise(across(tax_order, n_distinct))
  examining_dups6 <- examining_dups6 %>% rename(taxlineage6_n = tax_order)
  examining_dups7 <- thresholded_hit_nofish_dups %>% group_by(query) %>% summarise(across(tax_family, n_distinct))
  examining_dups7 <- examining_dups7 %>% rename(taxlineage7_n = tax_family)
  examining_dups8 <- thresholded_hit_nofish_dups %>% group_by(query) %>% summarise(across(tax_genus, n_distinct))
  examining_dups8 <- examining_dups8 %>% rename(taxlineage8_n = tax_genus)
  
  examining_dups_all <- mutate(examining_dups1, examining_dups2, examining_dups3, examining_dups4, examining_dups5, examining_dups6, examining_dups7, examining_dups8)
  rm(examining_dups1)
  rm(examining_dups2)
  rm(examining_dups3)
  rm(examining_dups4)
  rm(examining_dups5)
  rm(examining_dups6)
  rm(examining_dups7)
  rm(examining_dups8)
  
  thresholded_hit_nofish_dups2 <- full_join(thresholded_hit_nofish_dups,examining_dups_all)
  rm(thresholded_hit_nofish_dups)
  
  thresholded_hit_nofish_dups2$taxname_lca <- ifelse((thresholded_hit_nofish_dups2$taxlineage1_n > 1), "root", 
                                                     ifelse((thresholded_hit_nofish_dups2$taxlineage2_n > 1), thresholded_hit_nofish_dups2$tax_superkingdom, 
                                                            ifelse((thresholded_hit_nofish_dups2$taxlineage3_n > 1), thresholded_hit_nofish_dups2$tax_clade,
                                                                   ifelse((thresholded_hit_nofish_dups2$taxlineage4_n > 1), thresholded_hit_nofish_dups2$tax_kingdom,
                                                                          ifelse((thresholded_hit_nofish_dups2$taxlineage5_n > 1), thresholded_hit_nofish_dups2$tax_phylum,
                                                                                 ifelse((thresholded_hit_nofish_dups2$taxlineage6_n > 1), thresholded_hit_nofish_dups2$tax_class,
                                                                                        ifelse((thresholded_hit_nofish_dups2$taxlineage7_n > 1), thresholded_hit_nofish_dups2$tax_order,
                                                                                               ifelse((thresholded_hit_nofish_dups2$taxlineage8_n > 1), thresholded_hit_nofish_dups2$tax_family, thresholded_hit_nofish_dups2$tax_genus))))))))
  
  thresholded_hit_nofish_dups2 <- thresholded_hit_nofish_dups2 %>% relocate(taxname_lca, .after = taxname)
  

  thresholded_hit_nofish_dups2$taxname_lca <- ifelse((is.na(thresholded_hit_nofish_dups2$taxname_lca) == TRUE), thresholded_hit_nofish_dups2$tax_superkingdom, thresholded_hit_nofish_dups2$taxname_lca)

  ## extract just query and taxname new, then rename that...then add to thresholded_hit_nofish_rest to make new thresholded_hit_nofish
  thresholded_hit_nofish_dups3 <- thresholded_hit_nofish_dups2 %>% select(query,taxname_lca)
  thresholded_hit_nofish_dups4 <- thresholded_hit_nofishnomammals7 %>% group_by(query) %>% add_tally() %>% dplyr::filter(n > 1)
  thresholded_hit_nofish_dups3 <- thresholded_hit_nofish_dups3 %>% group_by(query) %>% slice(1)
  thresholded_hit_nofish_dups5 <- left_join(thresholded_hit_nofish_dups4, thresholded_hit_nofish_dups3) %>% 
    select(-taxname) %>% group_by(query) %>% slice(1) %>% relocate(taxname_lca, .after = taxid) %>% select(-taxoncategory)
  rm(thresholded_hit_nofish_dups3)
  rm(thresholded_hit_nofish_dups4)
  
  ## we need to refresh the _rest file!!
  thresholded_hit_nofish_rest <- thresholded_hit_nofish_rest %>% rename(taxname_lca = taxname) %>% select(-taxoncategory)
  
  
  ## MAY NEED TO DROP 'n' in first 2 cases??? select everything but n just in case
  if (nrow(thresholded_hit_nofish_rest) == 0) {
    thresholded_hit_nofishnomammals <- thresholded_hit_nofish_dups5 #%>% select(runname,query,node,length,coverage,rel_abundance,target,taxid,taxname_lca,bits,superkingdom,tax_kingdom,phylum,tax_class,tax_order,tax_family,tax_genus)
  } else if (nrow(thresholded_hit_nofish_dups5) == 0) {
    thresholded_hit_nofishnomammals <- thresholded_hit_nofish_rest #%>% select(runname,query,node,length,coverage,rel_abundance,target,taxid,taxname_lca,bits,superkingdom,tax_kingdom,phylum,tax_class,tax_order,tax_family,tax_genus)
  } else {
    thresholded_hit_nofishnomammals <- bind_rows(thresholded_hit_nofish_rest, thresholded_hit_nofish_dups5) %>% 
      select(-n) #%>% relocate(rel_abundance, .after = coverage)
  }
  
  ## add double-check to make sure n column is removed in edge cases above
  if ("n" %in% colnames(thresholded_hit_nofishnomammals)) {
    thresholded_hit_nofishnomammals <- thresholded_hit_nofishnomammals %>% select(-n)
  }
  
  ## TWEAKED TAXONOMIZR CODE SO WE CAN SEARCH FOR SAR 
  
  thresholded_hit_nofishnomammals$taxoncategory2 <- ifelse((grepl("root", thresholded_hit_nofishnomammals$taxname_lca) == TRUE), "Root_unresolved",
                                                           ifelse((grepl("Archaea", thresholded_hit_nofishnomammals$tax_superkingdom) == TRUE), "Archaea",
                                                                  ifelse((grepl("Viruses", thresholded_hit_nofishnomammals$tax_superkingdom) == TRUE), "Viruses",
                                                                         ifelse((grepl("Bacteria", thresholded_hit_nofishnomammals$tax_superkingdom) == TRUE), "Bacteria", 
                                                                                ifelse((grepl("Fungi", thresholded_hit_nofishnomammals$tax_kingdom) == TRUE), "Fungi", 
                                                                                       ifelse((grepl("Viridiplantae", thresholded_hit_nofishnomammals$tax_kingdom) == TRUE), "Plants", 
                                                                                              ifelse((grepl("Sar", thresholded_hit_nofishnomammals$tax_clade) == TRUE), "SAR_Eukaryotes",
                                                                                                     ifelse((grepl("Actinopteri", thresholded_hit_nofishnomammals$tax_class) == TRUE), "Fishes",
                                                                                                            ifelse((grepl("Primates", thresholded_hit_nofishnomammals$tax_order) == TRUE), "Human",
                                                                                                                   ifelse((thresholded_hit_nofishnomammals$taxname_lca == "Eukaryota"), "unresolved_Eukaryota", thresholded_hit_nofishnomammals$tax_phylum))))))))))
  
  ### CHANGING THIS AND BELOW NA FROM LCA TO "Root_unresolved" since these are usually artificial, plasmids, unidentified, etc.
  thresholded_hit_nofishnomammals <- thresholded_hit_nofishnomammals %>% mutate(across(taxoncategory2, ~ replace_na(., "Root_unresolved")))
  rm(thresholded_hit_nofish_dups5)
  rm(thresholded_hit_nofish_rest)
  rm(examining_dups_all)
  
  ## REMOVE TARGET_TITLE & ENSURE ORDER OF 4 COLUMNS
  thresholded_hit_nofishnomammals <- thresholded_hit_nofishnomammals %>% relocate(target, .after = rel_abundance) %>% relocate(bits, .after = target) %>% relocate(taxid, .after = bits) %>% relocate(taxname_lca, .after = taxid)
  if("target_title" %in% colnames(thresholded_hit_nofishnomammals) == FALSE) {
    thresholded_hit_nofishnomammals$target_title <- "N/A"
    thresholded_hit_nofishnomammals <- thresholded_hit_nofishnomammals %>% relocate(target_title, .before = taxoncategory2)
  }
  
  ## rename columns 8-36 by appending _NR
  thresholded_hit_nofishnomammals_diamond <- thresholded_hit_nofishnomammals %>% rename_with(~paste0(.x, "_NR"), c(target:taxoncategory2))
  
  currentchunkname <- chunknames[i]
  
  ## save then remove variables
  write.table(thresholded_hit_nofishnomammals_diamond, file = paste0(currentchunkname,"_diamond_hits_nofishnohuman",".tab"), sep = "\t", row.names = FALSE, quote = FALSE)
  rm(thresholded_hit_nofishnomammals)
  rm(thresholded_hit_nofishnomammals_diamond)
  
}
