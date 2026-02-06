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



## here would be a good place to run blastdbdmd command in slurm script then pick back up and run last parts as _lastbit.R script?
#allchunks_diamondnr_andblastntclustered <- read_tsv("taxonomy_hits_nonhost_mostrecent.tsv.gz")

#gzip -k taxonomy_hits_viruses0_fullcols_mostrecent.tsv
## then rename OG to _without_taxids and _with_taxids back to regular name
## then reload nonhost

allchunks_diamondnr_andblastntclustered <- read_tsv("taxonomy_hits_nonhost_mostrecent.tsv.gz")


## only save the undated version, then save the dated version in the withtaxid script!
write_tsv(allchunks_diamondnr_andblastntclustered, file = paste0("taxonomy_hits_nonhost_",Sys.Date(),".tsv.gz"))

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
write.table(allchunks_diamondnr_andblastntclustered_SAR_Eukaryotes, file = paste0("taxonomy_hits_SAReukaryotes_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_archaea, file = paste0("taxonomy_hits_archaea_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_mollusca, file = paste0("taxonomy_hits_mollusca_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_annelida, file = paste0("taxonomy_hits_annelida_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_nematoda, file = paste0("taxonomy_hits_nematoda_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_platyhelminthes, file = paste0("taxonomy_hits_platyhelminthes_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


## save lists for fasta pulling
#allchunks_diamondnr_andblastntclustered_viruses0 <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR == "Viruses")
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
write.table(allchunks_diamondnr_andblastntclustered_viruses_q, file = paste0("taxonomy_hits_viruses_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_bacteria_q, file = paste0("taxonomy_hits_bacteria_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_arthropoda_q, file = paste0("taxonomy_hits_arthropoda_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_plants_q, file = paste0("taxonomy_hits_plants_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_chordates_q, file = paste0("taxonomy_hits_chordates_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_fungi_q, file = paste0("taxonomy_hits_fungi_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_otherEukaryota_q, file = paste0("taxonomy_hits_otherEukaryota_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_SAR_Eukaryotes_q, file = paste0("taxonomy_hits_SAReukaryotes_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_archaea_q, file = paste0("taxonomy_hits_archaea_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_mollusca_q, file = paste0("taxonomy_hits_mollusca_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_annelida_q, file = paste0("taxonomy_hits_annelida_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_nematoda_q, file = paste0("taxonomy_hits_nematoda_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_platyhelminthes_q, file = paste0("taxonomy_hits_platyhelminthes_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

## also full list!
allchunks_diamondnr_andblastntclustered_q <- allchunks_diamondnr_andblastntclustered %>% select(query)
write.table(allchunks_diamondnr_andblastntclustered_q, file = paste0("taxonomy_hits_nonhost_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

## full list but without chordates
#allchunks_diamondnr_andblastntclustered <- read_tsv("taxonomy_hits_nonhost_mostrecent.tsv.gz")

allchunks_diamondnr_andblastntclustered_allbutchordates <- allchunks_diamondnr_andblastntclustered %>% dplyr::filter(taxoncategorysimple_NTorNR != "Chordata")

write_tsv(allchunks_diamondnr_andblastntclustered_allbutchordates, file = paste0("taxonomy_hits_nonhostnochordates_",Sys.Date(),".tsv.gz"))
write_tsv(allchunks_diamondnr_andblastntclustered_allbutchordates, file = paste0("taxonomy_hits_nonhostnochordates_mostrecent.tsv.gz"))
allchunks_diamondnr_andblastntclustered_allbutchordates_q <- allchunks_diamondnr_andblastntclustered_allbutchordates %>% select(query)
write.table(allchunks_diamondnr_andblastntclustered_allbutchordates_q, file = paste0("taxonomy_hits_nonhostnochordates_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


# ## seqtk commands will be in a separate script

# ####################################################################################
# ###### all curated virus code separated, will remove for pipeline
# ####################################################################################


# #######################################################################################
# #### outputing plots and table last

# ggsave(filename = paste("taxonomy_hits_nonhost_alluvialplot_all_",Sys.Date(),".png", sep=""), alluvial_plotall, width = 18, height = 9, units = "in", limitsize = FALSE)
# #ggsave(filename = paste("allchunks_alluvialplot_all0_",Sys.Date(),".png", sep=""), alluvial_plotall, width = 5.6, height = 3.4, units = "in", limitsize = FALSE)
# ggsave(filename = paste("taxonomy_hits_nonhost_alluvialplot_all_",Sys.Date(),".pdf", sep=""), alluvial_plotall, width = 18, height = 9, units = "in", limitsize = FALSE)

# write.table(data_alluviala, file = paste0("taxonomy_hits_nonhost_alluvialplot_counts_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# ggsave(filename = paste("taxonomy_hits_nonhost_treemap_",Sys.Date(),".png", sep=""), NTNRcontigs_treemap, width = 18, height = 9, units = "in", limitsize = FALSE)
# ggsave(filename = paste("taxonomy_hits_nonhost_treemap_",Sys.Date(),".pdf", sep=""), NTNRcontigs_treemap, width = 18, height = 9, units = "in", limitsize = FALSE)

# ggsave(filename = paste("taxonomy_hits_nonhost_treemap14categories_",Sys.Date(),".png", sep=""), NTNRcontigs_treemap14, width = 18, height = 9, units = "in", limitsize = FALSE)
# ggsave(filename = paste("taxonomy_hits_nonhost_treemap14categories_",Sys.Date(),".pdf", sep=""), NTNRcontigs_treemap14, width = 18, height = 9, units = "in", limitsize = FALSE)


# write.table(data_heatmap, file = paste0("taxonomy_hits_nonhost_treemap_counts_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)