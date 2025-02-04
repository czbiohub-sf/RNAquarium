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

write.table(thresholded_hit_nofishnomammalsNT, file = "allchunks_blastnclustered_hits_nonzfhum.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#######
## combine all NR outputs

#blastoutputsNR <- fs::dir_ls(diamondpath, glob="*diamond_hits_nonzfhum.tab", recurse = FALSE)
blastoutputsNR <- fs::dir_ls(diamondpath, glob="*.tab", recurse = FALSE)

#thresholded_hit_nofishnomammals <- read_tsv(blastoutputs)

thresholded_hit_nofishnomammalsNR <- read_tsv(blastoutputsNR)
head(thresholded_hit_nofishnomammalsNR)

write.table(thresholded_hit_nofishnomammalsNR, file = "allchunks_diamond_hits_nonzfhum.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

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


## also run a lowcoverage_flag
allchunks_diamondnr_andblastntclustered <-  allchunks_diamondnr_andblastntclustered %>%
  mutate(
    lowcoverage_flag = case_when(
      (sumperc_cov_NTclustered < 0.5 & sumperc_cov_NR < 0.5) ~ "NTNR_lowconfidence",
      sumperc_cov_NTclustered < 0.5 ~ "NT_lowconfidence",
      sumperc_cov_NR < 0.5 ~ "NR_lowconfidence",
      .default = ""
    )
  )

allchunks_diamondnr_andblastntclustered <- allchunks_diamondnr_andblastntclustered %>% relocate(lowcoverage_flag, .before = target_NTclustered)

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
    taxname_lca_NTorNR = if_else(!is.na(bits_NR) & !is.na(bits_NTclustered) & bits_NR > bits_NTclustered, taxname_lca_NR, taxname_lca_NTorNR),
    taxoncategory_NTorNR = if_else(!is.na(bits_NR) & !is.na(bits_NTclustered) & bits_NR > bits_NTclustered, taxoncategory_NR, taxoncategory_NTorNR),
    taxoncategorysimple_NTorNR = if_else(!is.na(bits_NR) & !is.na(bits_NTclustered) & bits_NR > bits_NTclustered, taxoncategorysimple_NR, taxoncategorysimple_NTorNR),
    bits_NTorNR = if_else(!is.na(bits_NR) & !is.na(bits_NTclustered) & bits_NR > bits_NTclustered, bits_NR, bits_NTorNR),
    evalue_NTorNR = if_else(!is.na(bits_NR) & !is.na(bits_NTclustered) & bits_NR > bits_NTclustered, evalue_NR, evalue_NTorNR)
  )

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
allchunks_diamondnr_andblastntclustered_NRbetter <- allchunks_diamondnr_andblastntclustered %>% filter(bits_NR > bits_NTclustered)
allchunks_diamondnr_andblastntclustered_taxonmismatch <- allchunks_diamondnr_andblastntclustered %>% filter(taxoncategorysimple_NR != taxoncategorysimple_NTclustered)
## then save
write.table(allchunks_diamondnr_andblastntclustered, file = "allchunks_blastnanddiamond_hits_nonzfhum.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_NRbetter, file = "allchunks_blastnanddiamond_hits_nonzfhum_NRbetter.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_taxonmismatch, file = "allchunks_blastnanddiamond_hits_nonzfhum_taxonmismatch.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## also always have a generic version saved NOTE THIS WILL CHANGE, BUT IS USED FOR PART 3 AND ALSO VIRUS CURATION...
write.table(allchunks_diamondnr_andblastntclustered, file = "allchunks_blastnanddiamond_hits_nonzfhum_mostrecent.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


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



## then save
write.table(allchunks_diamondnr_andblastntclustered_viruses, file = "allchunks_blastnanddiamond_hits_viruses.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## save all groups
write.table(allchunks_diamondnr_andblastntclustered_bacteria, file = "allchunks_blastnanddiamond_hits_bacteria.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_arthropoda, file = "allchunks_blastnanddiamond_hits_arthropoda.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_plants, file = "allchunks_blastnanddiamond_hits_plants.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_chordates, file = "allchunks_blastnanddiamond_hits_chordates.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_fungi, file = "allchunks_blastnanddiamond_hits_fungi.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_otherEukaryota, file = "allchunks_blastnanddiamond_hits_otherEukaryota.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_SAR_Eukaryotes, file = "allchunks_blastnanddiamond_hits_SAR_Eukaryotes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_archaea, file = "allchunks_blastnanddiamond_hits_archaea.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_mollusca, file = "allchunks_blastnanddiamond_hits_mollusca.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_annelida, file = "allchunks_blastnanddiamond_hits_annelida.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_nematoda, file = "allchunks_blastnanddiamond_hits_nematoda.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_platyhelminthes, file = "allchunks_blastnanddiamond_hits_platyhelminthes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


## save lists for fasta pulling
allchunks_diamondnr_andblastntclustered_viruses_q <- allchunks_diamondnr_andblastntclustered_viruses %>% select(query)
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
write.table(allchunks_diamondnr_andblastntclustered_viruses_q, file = "allchunks_blastnanddiamond_hits_viruses_list.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_bacteria_q, file = "allchunks_blastnanddiamond_hits_bacteria_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_arthropoda_q, file = "allchunks_blastnanddiamond_hits_arthropoda_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_plants_q, file = "allchunks_blastnanddiamond_hits_plants_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_chordates_q, file = "allchunks_blastnanddiamond_hits_chordates_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_fungi_q, file = "allchunks_blastnanddiamond_hits_fungi_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_otherEukaryota_q, file = "allchunks_blastnanddiamond_hits_otherEukaryota_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_SAR_Eukaryotes_q, file = "allchunks_blastnanddiamond_hits_SAR_Eukaryotes_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_archaea_q, file = "allchunks_blastnanddiamond_hits_archaea_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_mollusca_q, file = "allchunks_blastnanddiamond_hits_mollusca_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_annelida_q, file = "allchunks_blastnanddiamond_hits_annelida_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_nematoda_q, file = "allchunks_blastnanddiamond_hits_nematoda_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_platyhelminthes_q, file = "allchunks_blastnanddiamond_hits_platyhelminthes_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)


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
#fasta_viruses <- read.fasta("allchunks_blastnanddiamond_hits_viruses_list.fasta")
## rename seq.name to fullquery
#fasta_viruses <- fasta_viruses %>% rename(query = seq.name)

#allchunks_diamondnr_andblastntclustered_viruses_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_viruses, fasta_viruses)
#write.table(allchunks_diamondnr_andblastntclustered_viruses_withsequence, file = "RNaquarium_allchunks_blastnanddiamond_hits_viruses_withsequence.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



####################################################################################
###### alluvial code
####################################################################################

# thresholded_hit_nofishnomammalsNTsm <- thresholded_hit_nofishnomammalsNT %>% select(query,taxoncategory_NTclustered)
# thresholded_hit_nofishnomammalsNRsm <- thresholded_hit_nofishnomammalsNR %>% select(query,taxoncategory_NR)
# allchunks_diamondnr_andblastntclustered_lists <- full_join(thresholded_hit_nofishnomammalsNTsm,thresholded_hit_nofishnomammalsNRsm)

allchunks_diamondnr_andblastntclustered_lists <- allchunks_diamondnr_andblastntclustered %>% select(query,taxoncategorysimple_NTorNR,taxoncategorysimple_NTclustered,taxoncategorysimple_NR)


## numbers


## now new alluvial plots...
library(ggalluvial)

data_alluviala <- allchunks_diamondnr_andblastntclustered_lists %>%
  group_by(taxoncategorysimple_NTclustered, taxoncategorysimple_NR) %>%
  summarise(count = n()) %>%
  ungroup()

# Replace NA values with 'Missing'
data_alluviala <- data_alluviala %>%
  mutate(across(everything(), ~replace_na(.x, "Missing")))

## change theme_minimal to theme_classic or theme_gray
#alluvial_plotv0 <- alluvial_plotv data_alluvial
# alluvial_plotall <- ggplot(data_alluviala, aes(axis1 = taxoncategorysimple_NTclustered, axis2 = taxoncategorysimple_NR, y = count)) +
#   geom_alluvium(aes(fill = taxoncategorysimple_NTclustered), width = 1/2) +
#   geom_stratum() +
#   geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("taxoncategorysimple_NTclustered", "taxoncategorysimple_NR"), expand = c(0.08, 0.05)) +
#   theme_grey() +
#   labs(title = "Alluvial Diagram of NT vs NR taxonomic categories",
#        x = "",
#        y = "Count") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 

## stat_stratum(decreasing = TRUE) +
alluvial_plotall <- ggplot(data_alluviala, aes(axis1 = taxoncategorysimple_NTclustered, axis2 = taxoncategorysimple_NR, y = count)) +
  geom_alluvium(aes(fill = taxoncategorysimple_NTclustered), width = 1/2, decreasing = FALSE) +
  stat_stratum(decreasing = FALSE) +
  stat_stratum(geom = "text", aes(label = after_stat(stratum)), decreasing = FALSE, size = 2.5, min.y = 100) +
  scale_x_discrete(limits = c("taxoncategorysimple_NTclustered", "taxoncategorysimple_NR"), expand = c(0.08, 0.05)) +
#  scale_y_continuous(transform = "pseudo_log") +
#  ggfittext::geom_fit_text(stat = "stratum", width = 1/4, min.size = 2, label = "") +
  theme_grey() +
  labs(title = "Alluvial Diagram of NT vs NR taxonomic categories",
       x = "",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#alluvial_plotall

ggsave(filename = "allchunks_alluvialplot_all.png", alluvial_plotall, width = 17, height = 28, units = "in", limitsize = FALSE)
ggsave(filename = "allchunks_alluvialplot_all.pdf", alluvial_plotall, width = 17, height = 28, units = "in", limitsize = FALSE)


write.table(data_alluviala, file = "allchunks_alluvialplot_counts.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



####################################################################################
###### get numbers for heatmaps
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
data_heatmap2 <- data_heatmap %>% arrange(desc(count)) %>% slice_head(n = 10)


## new column combining text & numbers, also rewording unite & str_
data_heatmap2 <- data_heatmap2 %>% unite(label, taxoncategorysimple_NTorNR, count, sep = "\n", remove = FALSE)
#data_heatmap2 <- data_heatmap2 %>% unite(label, taxoncategorysimple_NTorNR, count, sep = ", ", remove = FALSE)
data_heatmap2 <- data_heatmap2 %>% mutate(label = str_replace_all(label, c("other_Eukaryota" = "other Eukaryota")))
data_heatmap2 <- data_heatmap2 %>% mutate(label = str_replace_all(label, c("Root_unresolved" = "Root (unresolved)")))

library(treemapify)
## 
# Plotting TreeMap Graph 
NTNRcontigs_treemap <- ggplot2::ggplot(data_heatmap2,aes(area=count,fill=taxoncategorysimple_NTorNR,label=label,subgroup=taxoncategorysimple_NTorNR)) + 
  treemapify::geom_treemap(layout="squarified") + 
  geom_treemap_text(place = "centre", size = 18) + 
  labs(title="Treemap of contigs found by NT + NR searches", fill="Taxon")


ggsave(filename = "RNaquarium_allchunks_blastnanddiamond_allhits_treemap.png", NTNRcontigs_treemap, width = 18, height = 9, units = "in", limitsize = FALSE)
ggsave(filename = "RNaquarium_allchunks_blastnanddiamond_allhits_treemap.pdf", NTNRcontigs_treemap, width = 18, height = 9, units = "in", limitsize = FALSE)




####################################################################################
###### all curated virus code separated, will remove for pipeline
####################################################################################
