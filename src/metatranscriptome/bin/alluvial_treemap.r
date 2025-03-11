#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(scales)
library(fs)
library(taxonomizr)
library(tidyverse)

parser <- ArgumentParser(description = "Alluvial Treemap Plotter")
parser$add_argument("-nt", "--nt-blast-dir")
parser$add_argument("-nr", "--nr-diamond-dir")
args <- parser$parse_args()

blastpath <- args$nt_blast_dir
diamondpath <- args$nr_diamond_dir

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


### combine all NT outputs, then subset viruses

#blastoutputs <- fs::dir_ls(blastpath, glob="blastnclustered_hits_nonzfhum", recurse = FALSE)
#virusoutputs <- fs::dir_ls(workingpath, glob="*tsv", recurse = FALSE)
#thresholded_hit_viruses_all <- read_tsv(virusoutputs)


## combining all NT outputs
#blastoutputsNT <- fs::dir_ls(blastpath, glob="*blastnclustered_hits_nonzfhum_2024-09-30.tab", recurse = FALSE)
blastoutputsNT <- fs::dir_ls(blastpath, glob="*.tab", recurse = FALSE)

thresholded_hit_nofishnomammalsNT <- read_tsv(blastoutputsNT)
head(thresholded_hit_nofishnomammalsNT)

write.table(thresholded_hit_nofishnomammalsNT, file = "allchunks_blastn_hits_nonhost.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#######
## combine all NR outputs

#blastoutputsNR <- fs::dir_ls(diamondpath, glob="*diamond_hits_nonzfhum.tab", recurse = FALSE)
blastoutputsNR <- fs::dir_ls(diamondpath, glob="*.tab", recurse = FALSE)

#thresholded_hit_nofishnomammals <- read_tsv(blastoutputs)

thresholded_hit_nofishnomammalsNR <- read_tsv(blastoutputsNR)
head(thresholded_hit_nofishnomammalsNR)

write.table(thresholded_hit_nofishnomammalsNR, file = "allchunks_diamond_hits_nonhost.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

####################################################################################














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






### then relocate
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxname_lca_NTclustered, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategory_NTclustered, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategorysimple_NTclustered, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxname_lca_NR, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategory_NR, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategorysimple_NTclustered, .before = target_NTclustered)
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(taxoncategorysimple_NR, .before = target_NTclustered)


### then coalesces here

## update - we will want to make rules so that if NT bitscore is much lower than NR, to use NR for all fields
### also note commands for BLAST step for virus counts are below search DECEMBER UPDATE

### rule would go like this:
# first the coalesce as before to replace all missing NT with NR
## THEN
## only when bits_NTclustered < bits_NR replace those 4 created NTorNR fields with the NR field instead...

# allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% rename(taxoncategory_NR = taxoncategory2_NR)
# allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% rename(taxoncategory_NTclustered = taxoncategory2_NTclustered)
# allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% rename(taxoncategory_NTorNR = taxoncategory2_NTorNR)
# allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% rename(taxoncategorysimple_NR = taxoncategory2simple_NR)
# allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% rename(taxoncategorysimple_NTclustered = taxoncategory2simple_NTclustered)
# allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% rename(taxoncategorysimple_NTorNR = taxoncategory2simple_NTorNR)

## rename

allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% mutate(taxname_lca_NTorNR = coalesce(taxname_lca_NTclustered,taxname_lca_NR))
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% mutate(taxoncategory_NTorNR = coalesce(taxoncategory_NTclustered,taxoncategory_NR))
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% mutate(taxoncategorysimple_NTorNR = coalesce(taxoncategorysimple_NTclustered,taxoncategorysimple_NR))


allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% mutate(bits_NTorNR = coalesce(bits_NTclustered,bits_NR))
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% mutate(evalue_NTorNR = coalesce(evalue_NTclustered,evalue_NR))



allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>%
  mutate(
    taxname_lca_NTorNR = if_else(!is.na(bits_NR) & !is.na(bits_NTclustered) & bits_NR >= bits_NTclustered, taxname_lca_NR, taxname_lca_NTorNR),
    taxoncategory_NTorNR = if_else(!is.na(bits_NR) & !is.na(bits_NTclustered) & bits_NR >= bits_NTclustered, taxoncategory_NR, taxoncategory_NTorNR),
    taxoncategorysimple_NTorNR = if_else(!is.na(bits_NR) & !is.na(bits_NTclustered) & bits_NR >= bits_NTclustered, taxoncategorysimple_NR, taxoncategorysimple_NTorNR),
    bits_NTorNR = if_else(!is.na(bits_NR) & !is.na(bits_NTclustered) & bits_NR >= bits_NTclustered, bits_NR, bits_NTorNR),
    evalue_NTorNR = if_else(!is.na(bits_NR) & !is.na(bits_NTclustered) & bits_NR >= bits_NTclustered, evalue_NR, evalue_NTorNR)
  )

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

allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(analysis_used, .before = taxname_lca_NTorNR)



## finally add the lowcoverage flag only on the analysis_used
## also run a lowcoverage_flag
allchunks_diamondnr_andblastntclustered <-  allchunks_diamondnr_andblastntclustered %>%
  mutate(
    lowcoverage_flag = case_when(
      (analysis_used == "NR" & sumperc_cov_NR < 0.5) ~ "NR_lowconfidence",
      (analysis_used == "NTclustered" & sumperc_cov_NTclustered < 0.5) ~ "NT_lowconfidence",
      .default = ""
    )
  )

allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(lowcoverage_flag, .before = target_NTclustered)


## relocates last
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

# get mismatches...
allchunks_diamondnr_andblastntclustered_NRbetter <- allchunks_diamondnr_andblastntclustered %>% filter(bits_NR >= bits_NTclustered)
allchunks_diamondnr_andblastntclustered_taxonmismatch <- allchunks_diamondnr_andblastntclustered %>% filter(taxoncategorysimple_NR != taxoncategorysimple_NTclustered)
## then save - NOTE BEFORE REMOVING COLUMNS WE ARE SAVING ONLY THESE 3

#write.table(allchunks_diamondnr_andblastntclustered, file = paste0("allchunks_blastnanddiamond_hits_nofishnohuman_fullcols_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(allchunks_diamondnr_andblastntclustered_NRbetter, file = paste0("allchunks_blastnanddiamond_hits_nofishnohuman_NRbetter_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(allchunks_diamondnr_andblastntclustered_taxonmismatch, file = paste0("allchunks_blastnanddiamond_hits_nofishnohuman_taxonmismatch_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


## then save
write.table(allchunks_diamondnr_andblastntclustered, file = "taxonomy_hits_nonhost_fullcols_mostrecent.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_NRbetter, file = "taxonomy_hits_nonhost_NRbetter.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_taxonmismatch, file = "taxonomy_hits_nonhost_taxonmismatch.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## also reverse select fields based on analysis_used! do this at very end after alluvial plots & treemaps
## but save viruses0 file before removing these fields...they are used for virus steps
## also save a full_cols version
## now save ONLY viruses version - viruses00 is only used for alluvial plots
allchunks_diamondnr_andblastntclustered_viruses0 <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "Viruses")
write.table(allchunks_diamondnr_andblastntclustered_viruses0, file = "taxonomy_hits_viruses_fullcols_mostrecent.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

allchunks_diamondnr_andblastntclustered_viruses00 <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategory_NTclustered == "Viruses" | taxoncategory_NR == "Viruses")
write.table(allchunks_diamondnr_andblastntclustered_viruses00, file = "taxonomy_hits_viruses0_fullcols_mostrecent.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


## now alluvial plots & treemaps
####################################################################################
###### alluvial code
####################################################################################

# thresholded_hit_nofishnomammalsNTsm <- thresholded_hit_nofishnomammalsNT %>% select(query,taxoncategory_NTclustered)
# thresholded_hit_nofishnomammalsNRsm <- thresholded_hit_nofishnomammalsNR %>% select(query,taxoncategory_NR)
# allchunks_diamondnr_andblastntclustered_lists <- full_join(thresholded_hit_nofishnomammalsNTsm,thresholded_hit_nofishnomammalsNRsm)
##write.table(allchunks_diamondnr_andblastntclustered, file = paste0("taxonomy_hits_nonhost_fullcols_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
#allchunks_diamondnr_andblastntclustered <- read_tsv("taxonomy_hits_nonhost_fullcols_mostrecent.tsv")
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
#data_alluviala <- read_tsv("taxonomy_hits_nonhost_alluvialplot_counts.tsv")
data_alluviala <- data_alluviala %>%
  mutate(across(everything(), ~replace_na(.x, "Missing")))

### moving saving of plots and table to very end
#write.table(data_alluviala, file = paste0("taxonomy_hits_nonhost_alluvialplot_counts_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

data_alluviala <- data_alluviala %>% rename(NT = taxoncategorysimple_NTclustered)
data_alluviala <- data_alluviala %>% rename(NR = taxoncategorysimple_NR)
data_alluviala <- data_alluviala %>% rename(used_NTorNR = taxoncategorysimple_NTorNR)

sumcountonepercent <- (sum(data_alluviala$count) / 100)

## set of order of each category using a factor??
data_alluviala <- data_alluviala %>% arrange(desc(count))
alluvialorder2 <- unique(data_alluviala$NT)


## remove rows below cutoff of sumcountonepercent
sumcountonepercent4 <- sumcountonepercent / 4
data_alluvial <- data_alluviala %>% dplyr::filter(count > sumcountonepercent4)


## edge case when there are groups not in NT, need to combine with other columns
##data_alluviala %>% expand(NT, used_NTorNR, NR) use complete - instead use pivot_longer
alluvialorder <- data_alluvial %>% pivot_longer(cols = -count)
alluvialorder <- alluvialorder %>% arrange(desc(name),desc(count))
alluvialorder2 <- unique(alluvialorder$value)
## make sure there are only 9 categories, per Set1 paletted??
#alluvialorder2 <- head(alluvialorder2, 9)

data_alluvial$NT <- factor(data_alluvial$NT, levels = alluvialorder2, ordered = TRUE)
data_alluvial$NR <- factor(data_alluvial$NR, levels = alluvialorder2, ordered = TRUE)
data_alluvial$used_NTorNR <- factor(data_alluvial$used_NTorNR, levels = alluvialorder2, ordered = TRUE)
## this works, but also have to remove "decreasing = FALSE" across the commands


## stat_stratum(decreasing = TRUE) +
# alluvial_plotall <- ggplot(data_alluvial, aes(axis1 = NT, axis2 = NR, y = count)) +
#   geom_alluvium(aes(fill = NT), width = 1/2, decreasing = FALSE) +
#   stat_stratum(decreasing = FALSE) +
#   stat_stratum(geom = "text", aes(label = after_stat(stratum)), decreasing = FALSE, size = 5, min.y = sumcountonepercent) +
#   scale_x_discrete(limits = c("NT", "NR"), expand = c(0.08, 0.05)) +
#   #  scale_y_continuous(transform = "pseudo_log") +
#   #  ggfittext::geom_fit_text(stat = "stratum", width = 1/4, min.size = 2, label = "") +
#   theme_grey(base_family="Helvetica", base_size = 16) +
#   labs(title = "Alluvial Diagram of NT vs NR taxonomic categories",
#        x = "",
#        y = "Count") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face="bold"))

# data_alluvial <- data_alluvial %>%
#   mutate(
#     formatted_count = paste0((comma(count)))
#   )



alluvial_plotall <- ggplot(data_alluvial, aes(axis1 = NT, axis2 = used_NTorNR, axis3 = NR, y = count)) +
  geom_alluvium(aes(fill = used_NTorNR), width = 1/2) +
  stat_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 5, min.y = sumcountonepercent) +
  scale_x_discrete(limits = c("NT", "best (NT or NR)", "NR"), expand = c(0.08, 0.02)) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_brewer(palette = "Set1") + ## Set3 has more categories but worse color scheme
  theme_minimal(base_family="Helvetica", base_size = 16) + ## also theme_minimal, theme_void
#  labs(title = "Alluvial Diagram of NT vs NR taxonomic categories", fill = "best (NT or NR)",
  labs(title = "Alluvial Diagram of NT vs NR taxonomic categories",
       x = "",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16))

#alluvial_plotall
#alluvial_plotall , base_size = 16 family base_family="Helvetica", base_size = 16


### moving saving of plots to very end
# ggsave(filename = paste("taxonomy_hits_nonhost_alluvialplot_all_",Sys.Date(),".png", sep=""), alluvial_plotall, width = 18, height = 9, units = "in", limitsize = FALSE)
# #ggsave(filename = paste("allchunks_alluvialplot_all0_",Sys.Date(),".png", sep=""), alluvial_plotall, width = 5.6, height = 3.4, units = "in", limitsize = FALSE)
# ggsave(filename = paste("taxonomy_hits_nonhost_alluvialplot_all_",Sys.Date(),".pdf", sep=""), alluvial_plotall, width = 18, height = 9, units = "in", limitsize = FALSE)

#write.table(data_alluviala, file = paste0("taxonomy_hits_nonhost_alluvialplot_counts_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)



####################################################################################
###### get numbers for treemaps
####################################################################################


## coalesce bits for NT & NR
#allchunks_diamondnr_andblastntclustered_lists <- allchunks_diamondnr_andblastntclustered_lists %>% mutate(taxoncategorysimple_NTorNR = coalesce(taxoncategorysimple_NTclustered,taxoncategorysimple_NR))
# now already done for taxoncategorysimple_NTorNR
## then get numbers for each of the unique values of taxoncategorysimple_NTorNR
### group by then count
data_heatmap <- allchunks_diamondnr_andblastntclustered_lists %>%
  group_by(taxoncategorysimple_NTorNR) %>%
  summarise(count = n()) %>%
  ungroup()

## then sort descending, take just first 10, and use these for treemap
data_heatmap <- data_heatmap %>% arrange(desc(count))
data_heatmap2 <- data_heatmap %>% arrange(desc(count)) %>% slice_head(n = 10)

## save these numbers (now at end)
#write.table(data_heatmap, file = paste0("taxonomy_hits_nonhost_treemap_counts_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


## new column combining text & numbers, also rewording unite & str_
data_heatmap2 <- data_heatmap2 %>% unite(label, taxoncategorysimple_NTorNR, count, sep = "\n", remove = FALSE)
#data_heatmap2 <- data_heatmap2 %>% unite(label, taxoncategorysimple_NTorNR, count, sep = ", ", remove = FALSE)
data_heatmap2 <- data_heatmap2 %>% mutate(label = str_replace_all(label, c("other_Eukaryota" = "other Eukaryota")))
data_heatmap2 <- data_heatmap2 %>% mutate(label = str_replace_all(label, c("Root_unresolved" = "Root (unresolved)")))

library(treemapify)
## 
# Plotting TreeMap Graph 
# NTNRcontigs_treemap <- ggplot2::ggplot(data_heatmap2,aes(area=count,fill=taxoncategorysimple_NTorNR,label=label,subgroup=taxoncategorysimple_NTorNR)) + 
#   treemapify::geom_treemap(layout="squarified") + 
#   geom_treemap_text(place = "centre", size = 18, fontface = "italic") + 
#   labs(title="Treemap of contigs found by NT + NR searches", fill="Taxon")

## updating to get commas in numbers & italics
data_heatmap2 <- data_heatmap2 %>%
  mutate(
    formatted_label = paste0(taxoncategorysimple_NTorNR, "\n", (comma(count)))
  )

NTNRcontigs_treemap <- ggplot2::ggplot(data_heatmap2,aes(area=count,fill=taxoncategorysimple_NTorNR,label=formatted_label,subgroup=taxoncategorysimple_NTorNR)) + 
  treemapify::geom_treemap(layout="squarified") + 
  geom_treemap_text(place = "centre", size = 18, fontface = "italic") + 
  labs(title="Treemap of contigs found by NT + NR searches", fill="Taxon")

### moving saving of plots and table to very end
# ggsave(filename = paste("taxonomy_hits_nonhost_treemap_",Sys.Date(),".png", sep=""), NTNRcontigs_treemap, width = 18, height = 9, units = "in", limitsize = FALSE)
# ggsave(filename = paste("taxonomy_hits_nonhost_treemap_",Sys.Date(),".pdf", sep=""), NTNRcontigs_treemap, width = 18, height = 9, units = "in", limitsize = FALSE)


####################################################################################
####################################################################################

## now reverse select fields based on analysis_used! do this at very end after alluvial plots & treemaps
## not reverse select, just gsub those other fields to "" - can skip since we consolidate then remove...
# allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>%
#   mutate(
#     across(
#       matches("_NTclustered$"),
#       ~ if_else(analysis_used == "NR", NA, .)
#     ),
#     across(
#       matches("_NR$"),
#       ~ if_else(analysis_used == "NTclustered", NA, .)
#     )
#   )

## then also consolidate many columns, then just remove a bunch

allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>%
  mutate(
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

#allchunks_diamondnr_andblastntclustered_viruses <- allchunks_diamondnr_andblastntclustered_viruses %>% relocate(tax_superkingdom_NTorNR, .after = taxname_lca_NTorNR)

## then negative select
allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% select(-target_NR) %>% select(-taxid_NTclustered) %>% select(-taxid_NR) %>%
  select(-gene_NTclustered) %>% select(-gene_NR) %>% select(-allele_NTclustered) %>% select(-allele_NR) %>% select(-pident_NTclustered) %>% select(-pident_NR) %>% select(-sumperc_cov_NTclustered) %>% select(-sumperc_cov_NR) %>%
  select(-alnlen_NTclustered) %>% select(-alnlen_NR) %>% select(-mismatch_NTclustered) %>% select(-mismatch_NR) %>% select(-qcov_NTclustered) %>% select(-qcov_NR) %>% select(-gapopen_NTclustered) %>% select(-gapopen_NR) %>%
  select(-qstart_NTclustered) %>% select(-qstart_NR) %>% select(-qend_NTclustered) %>% select(-qend_NR) %>% select(-tstart_NTclustered) %>% select(-tstart_NR) %>% select(-tend_NTclustered) %>% select(-tend_NR) %>% 
  select(-target_title_NTclustered) %>% select(-target_title_NR) %>% select(-analysis_NTclustered) %>% select(-analysis_NR) %>% select(-sumalnlen_NTclustered) %>% select(-sumalnlen_NR) %>% select(-maxbits_NTclustered) %>% select(-maxbits_NR) %>%
  select(-bits_percmax_NTclustered) %>% select(-bits_percmax_NR) %>% select(-tax_superkingdom_NTclustered) %>% select(-tax_superkingdom_NR) %>% select(-tax_clade_NTclustered) %>% select(-tax_clade_NR) %>% select(-tax_kingdom_NTclustered) %>% select(-tax_kingdom_NR) %>%
  select(-tax_phylum_NTclustered) %>% select(-tax_phylum_NR) %>% select(-tax_class_NTclustered) %>% select(-tax_class_NR) %>% select(-tax_order_NTclustered) %>% select(-tax_order_NR) %>% select(-tax_family_NTclustered) %>% select(-tax_family_NR) %>%
  select(-tax_genus_NTclustered) %>% select(-tax_genus_NR) %>% select(-tax_species_NTclustered) %>% select(-tax_species_NR)

####################################################################################
####################################################################################

## now saving & subsetting
## then save
write.table(allchunks_diamondnr_andblastntclustered, file = paste0("taxonomy_hits_nonhost_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## also always have a generic version saved NOTE THIS WILL CHANGE, BUT IS USED FOR PART 3 AND ALSO VIRUS CURATION...
write.table(allchunks_diamondnr_andblastntclustered, file = paste0("taxonomy_hits_nonhost_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


## also ALL BROAD CATEGORIES ## revise code throughout second half to use updated NTorNR categories instead of either NT or NR...
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

## then save - virus0 here because more columns will be made in virus-specific scripts downstream
write.table(allchunks_diamondnr_andblastntclustered_viruses, file = paste0("taxonomy_hits_viruses0_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## save all groups
write.table(allchunks_diamondnr_andblastntclustered_bacteria, file = paste0("taxonomy_hits_bacteria_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_arthropoda, file = paste0("taxonomy_hits_arthropoda_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_plants, file = paste0("taxonomy_hits_plants_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_chordates, file = paste0("taxonomy_hits_chordates_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_fungi, file = paste0("taxonomy_hits_fungi_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_otherEukaryota, file = paste0("taxonomy_hits_otherEukaryota_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_SAR_Eukaryotes, file = paste0("taxonomy_hits_SAR_Eukaryotes_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_archaea, file = paste0("taxonomy_hits_archaea_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_mollusca, file = paste0("taxonomy_hits_mollusca_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_annelida, file = paste0("taxonomy_hits_annelida_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_nematoda, file = paste0("taxonomy_hits_nematoda_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_platyhelminthes, file = paste0("taxonomy_hits_platyhelminthes_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


## save lists for fasta pulling
allchunks_diamondnr_andblastntclustered_viruses_q <- allchunks_diamondnr_andblastntclustered_viruses0 %>% select(query)
allchunks_diamondnr_andblastntclustered_bacteria_q <- allchunks_diamondnr_andblastntclustered_bacteria %>% select(query)
allchunks_diamondnr_andblastntclustered_arthropoda_q <- allchunks_diamondnr_andblastntclustered_arthropoda %>% select(query)
allchunks_diamondnr_andblastntclustered_plants_q <- allchunks_diamondnr_andblastntclustered_plants %>% select(query)
allchunks_diamondnr_andblastntclustered_chordates_q <- allchunks_diamondnr_andblastntclustered_chordates %>% select(query)
allchunks_diamondnr_andblastntclustered_fungi_q <- allchunks_diamondnr_andblastntclustered_fungi %>% select(query)
allchunks_diamondnr_andblastntclustered_otherEukaryota_q <- allchunks_diamondnr_andblastntclustered_otherEukaryota %>% select(query)
allchunks_diamondnr_andblastntclustered_SAR_Eukaryotes_q <- allchunks_diamondnr_andblastntclustered_SAR_Eukaryotes %>% select(query)
allchunks_diamondnr_andblastntclustered_archaea_q <- allchunks_diamondnr_andblastntclustered_archaea %>% select(query)
allchunks_diamondnr_andblastntclustered_mollusca_q <- allchunks_diamondnr_andblastntclustered_mollusca %>% select(query)
allchunks_diamondnr_andblastntclustered_annelida_q <- allchunks_diamondnr_andblastntclustered_annelida %>% select(query)
allchunks_diamondnr_andblastntclustered_nematoda_q <- allchunks_diamondnr_andblastntclustered_nematoda %>% select(query)
allchunks_diamondnr_andblastntclustered_platyhelminthes_q <- allchunks_diamondnr_andblastntclustered_platyhelminthes %>% select(query)

## then save as text file, for seqtk command
write.table(allchunks_diamondnr_andblastntclustered_viruses_q, file = paste0("taxonomy_hits_viruses_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_bacteria_q, file = paste0("taxonomy_hits_bacteria_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_arthropoda_q, file = paste0("taxonomy_hits_arthropoda_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_plants_q, file = paste0("taxonomy_hits_plants_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_chordates_q, file = paste0("taxonomy_hits_chordates_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_fungi_q, file = paste0("taxonomy_hits_fungi_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_otherEukaryota_q, file = paste0("taxonomy_hits_otherEukaryota_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_SAR_Eukaryotes_q, file = paste0("taxonomy_hits_SAR_Eukaryotes_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_archaea_q, file = paste0("taxonomy_hits_archaea_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_mollusca_q, file = paste0("taxonomy_hits_mollusca_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_annelida_q, file = paste0("taxonomy_hits_annelida_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_nematoda_q, file = paste0("taxonomy_hits_nematoda_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_platyhelminthes_q, file = paste0("taxonomy_hits_platyhelminthes_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

## also full list!
allchunks_diamondnr_andblastntclustered_q <- allchunks_diamondnr_andblastntclustered %>% select(query)
write.table(allchunks_diamondnr_andblastntclustered_q, file = paste0("taxonomy_hits_nonhost_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

## seqtk commands will be in a separate script
### seqtk script is here
# 
# #!/bin/bash
# module purge
# module load anaconda/2023.03
# conda activate seqtk
# 
#for fasta in ./paired_end/*.fasta; do
#  date
#  seqtk subseq "$fasta" allchunks_blastnanddiamond_hits_viruses_list.txt >> allchunks_blastnanddiamond_hits_viruses_list.fasta
#  seqtk subseq "$fasta" allchunks_blastnanddiamond_hits_bacteria_list.txt >> allchunks_blastnanddiamond_hits_bacteria_list.fasta
#  seqtk subseq "$fasta" allchunks_blastnanddiamond_hits_arthropoda_list.txt >> allchunks_blastnanddiamond_hits_arthropoda_list.fasta
#  seqtk subseq "$fasta" allchunks_blastnanddiamond_hits_plants_list.txt >> allchunks_blastnanddiamond_hits_plants_list.fasta
#  seqtk subseq "$fasta" allchunks_blastnanddiamond_hits_chordates_list.txt >> allchunks_blastnanddiamond_hits_chordates_list.fasta
#  seqtk subseq "$fasta" allchunks_blastnanddiamond_hits_fungi_list.txt >> allchunks_blastnanddiamond_hits_fungi_list.fasta
#  seqtk subseq "$fasta" allchunks_blastnanddiamond_hits_otherEukaryota_list.txt >> allchunks_blastnanddiamond_hits_otherEukaryota_list.fasta
#  seqtk subseq "$fasta" allchunks_blastnanddiamond_hits_SAR_Eukaryotes_list.txt >> allchunks_blastnanddiamond_hits_SAR_Eukaryotes_list.fasta
#  seqtk subseq "$fasta" allchunks_blastnanddiamond_hits_archaea_list.txt >> allchunks_blastnanddiamond_hits_archaea_list.fasta
#  seqtk subseq "$fasta" allchunks_blastnanddiamond_hits_mollusca_list.txt >> allchunks_blastnanddiamond_hits_mollusca_list.fasta
#  seqtk subseq "$fasta" allchunks_blastnanddiamond_hits_annelida_list.txt >> allchunks_blastnanddiamond_hits_annelida_list.fasta
#  seqtk subseq "$fasta" allchunks_blastnanddiamond_hits_nematoda_list.txt >> allchunks_blastnanddiamond_hits_nematoda_list.fasta
#  seqtk subseq "$fasta" allchunks_blastnanddiamond_hits_platyhelminthes_list.txt >> allchunks_blastnanddiamond_hits_platyhelminthes_list.fasta
#done

#for fasta in *.fasta; do
#  seqtk subseq "$fasta" allchunks_blastnclustered_hits_viruses_list.txt >> allchunks_blastnclustered_hits_viruses_list.fasta
#done


## eventually will want a separate script to combine fastq sequences with above files
## then also - scripts for mmseqs cluster + adding these cluster name +info to all files? also minimap2 for viewing


## MAYBE SAVE UP TO SEQTK COMMANDS AS .RDATA, THEN LOAD BACK FOR THESE COMMANDS??
## THESE ARE NOW IN THIRD SCRIPT...
#library(phylotools)
#fasta_viruses <- read.fasta("taxonomy_hits_viruses_list.fasta")
## rename seq.name to fullquery
#fasta_viruses <- fasta_viruses %>% rename(query = seq.name)

#allchunks_diamondnr_andblastntclustered_viruses_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_viruses, fasta_viruses)
#write.table(allchunks_diamondnr_andblastntclustered_viruses_withsequence, file = paste0("RNaquarium_allchunks_blastnanddiamond_hits_viruses_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

####################################################################################
###### all curated virus code separated, will remove for pipeline
####################################################################################


#######################################################################################
#### outputing plots and table last

ggsave(filename = paste("taxonomy_hits_nonhost_alluvialplot_all_",Sys.Date(),".png", sep=""), alluvial_plotall, width = 18, height = 9, units = "in", limitsize = FALSE)
#ggsave(filename = paste("allchunks_alluvialplot_all0_",Sys.Date(),".png", sep=""), alluvial_plotall, width = 5.6, height = 3.4, units = "in", limitsize = FALSE)
ggsave(filename = paste("taxonomy_hits_nonhost_alluvialplot_all_",Sys.Date(),".pdf", sep=""), alluvial_plotall, width = 18, height = 9, units = "in", limitsize = FALSE)

write.table(data_alluviala, file = paste0("taxonomy_hits_nonhost_alluvialplot_counts_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

ggsave(filename = paste("taxonomy_hits_nonhost_treemap_",Sys.Date(),".png", sep=""), NTNRcontigs_treemap, width = 18, height = 9, units = "in", limitsize = FALSE)
ggsave(filename = paste("taxonomy_hits_nonhost_treemap_",Sys.Date(),".pdf", sep=""), NTNRcontigs_treemap, width = 18, height = 9, units = "in", limitsize = FALSE)

write.table(data_heatmap, file = paste0("taxonomy_hits_nonhost_treemap_counts_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
