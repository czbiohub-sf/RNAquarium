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
#allchunks_diamondnr_andblastntclustered <- read_tsv("taxonomy_hits_nonhost_mostrecent.tsv")
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



## seqtk commands will be in a separate script - part 2


## now part 3 commands
library(phylotools)


## now combining each of the 13 nonhost subsets
#allchunks_blastnanddiamond_hits_viruses_list.fasta
#allchunks_blastnanddiamond_hits_fungi_list.fasta
#allchunks_blastnanddiamond_hits_archaea_list.fasta
#allchunks_blastnanddiamond_hits_platyhelminthes_list.fasta
#allchunks_blastnanddiamond_hits_mollusca_list.fasta
#allchunks_blastnanddiamond_hits_arthropoda_list.fasta
#allchunks_blastnanddiamond_hits_nematoda_list.fasta
#allchunks_blastnanddiamond_hits_annelida_list.fasta
#allchunks_blastnanddiamond_hits_chordates_list.fasta
#allchunks_blastnanddiamond_hits_bacteria_list.fasta
#allchunks_blastnanddiamond_hits_plants_list.fasta
#allchunks_blastnanddiamond_hits_otherEukaryota_list.fasta
#allchunks_blastnanddiamond_hits_SAR_Eukaryotes_list.fasta

fasta_viruses <- read.fasta("taxonomy_hits_viruses_list.fasta")
## rename seq.name to fullquery
fasta_viruses <- fasta_viruses %>% rename(query = seq.name)
allchunks_diamondnr_andblastntclustered_viruses_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_viruses, fasta_viruses)
write.table(allchunks_diamondnr_andblastntclustered_viruses_withsequence, file = paste0("taxonomy_hits_viruses0_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

fasta_fungi <- read.fasta("taxonomy_hits_fungi_list.fasta")
## rename seq.name to fullquery
fasta_fungi <- fasta_fungi %>% rename(query = seq.name)
allchunks_diamondnr_andblastntclustered_fungi_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_fungi, fasta_fungi)
write.table(allchunks_diamondnr_andblastntclustered_fungi_withsequence, file = paste0("taxonomy_hits_fungi_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

fasta_archaea <- read.fasta("taxonomy_hits_archaea_list.fasta")
## rename seq.name to fullquery
fasta_archaea <- fasta_archaea %>% rename(query = seq.name)
allchunks_diamondnr_andblastntclustered_archaea_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_archaea, fasta_archaea)
write.table(allchunks_diamondnr_andblastntclustered_archaea_withsequence, file = paste0("taxonomy_hits_archaea_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

fasta_platyhelminthes <- read.fasta("taxonomy_hits_platyhelminthes_list.fasta")
## rename seq.name to fullquery
fasta_platyhelminthes <- fasta_platyhelminthes %>% rename(query = seq.name)
allchunks_diamondnr_andblastntclustered_platyhelminthes_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_platyhelminthes, fasta_platyhelminthes)
write.table(allchunks_diamondnr_andblastntclustered_platyhelminthes_withsequence, file = paste0("taxonomy_hits_platyhelminthes_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

fasta_mollusca <- read.fasta("taxonomy_hits_mollusca_list.fasta")
## rename seq.name to fullquery
fasta_mollusca <- fasta_mollusca %>% rename(query = seq.name)
allchunks_diamondnr_andblastntclustered_mollusca_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_mollusca, fasta_mollusca)
write.table(allchunks_diamondnr_andblastntclustered_mollusca_withsequence, file = paste0("taxonomy_hits_mollusca_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

fasta_arthropoda <- read.fasta("taxonomy_hits_arthropoda_list.fasta")
## rename seq.name to fullquery
fasta_arthropoda <- fasta_arthropoda %>% rename(query = seq.name)
allchunks_diamondnr_andblastntclustered_arthropoda_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_arthropoda, fasta_arthropoda)
write.table(allchunks_diamondnr_andblastntclustered_arthropoda_withsequence, file = paste0("taxonomy_hits_arthropoda_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

fasta_nematoda <- read.fasta("taxonomy_hits_nematoda_list.fasta")
## rename seq.name to fullquery
fasta_nematoda <- fasta_nematoda %>% rename(query = seq.name)
allchunks_diamondnr_andblastntclustered_nematoda_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_nematoda, fasta_nematoda)
write.table(allchunks_diamondnr_andblastntclustered_nematoda_withsequence, file = paste0("taxonomy_hits_nematoda_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

fasta_annelida <- read.fasta("taxonomy_hits_annelida_list.fasta")
## rename seq.name to fullquery
fasta_annelida <- fasta_annelida %>% rename(query = seq.name)
allchunks_diamondnr_andblastntclustered_annelida_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_annelida, fasta_annelida)
write.table(allchunks_diamondnr_andblastntclustered_annelida_withsequence, file = paste0("taxonomy_hits_annelida_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

fasta_chordates <- read.fasta("taxonomy_hits_chordates_list.fasta")
## rename seq.name to fullquery
fasta_chordates <- fasta_chordates %>% rename(query = seq.name)
allchunks_diamondnr_andblastntclustered_chordates_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_chordates, fasta_chordates)
write.table(allchunks_diamondnr_andblastntclustered_chordates_withsequence, file = paste0("taxonomy_hits_chordates_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

fasta_bacteria <- read.fasta("taxonomy_hits_bacteria_list.fasta")
## rename seq.name to fullquery
fasta_bacteria <- fasta_bacteria %>% rename(query = seq.name)
allchunks_diamondnr_andblastntclustered_bacteria_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_bacteria, fasta_bacteria)
write.table(allchunks_diamondnr_andblastntclustered_bacteria_withsequence, file = paste0("taxonomy_hits_bacteria_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

fasta_plants <- read.fasta("taxonomy_hits_plants_list.fasta")
## rename seq.name to fullquery
fasta_plants <- fasta_plants %>% rename(query = seq.name)
allchunks_diamondnr_andblastntclustered_plants_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_plants, fasta_plants)
write.table(allchunks_diamondnr_andblastntclustered_plants_withsequence, file = paste0("taxonomy_hits_plants_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

fasta_otherEukaryota <- read.fasta("taxonomy_hits_otherEukaryota_list.fasta")
## rename seq.name to fullquery
fasta_otherEukaryota <- fasta_otherEukaryota %>% rename(query = seq.name)
allchunks_diamondnr_andblastntclustered_otherEukaryota_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_otherEukaryota, fasta_otherEukaryota)
write.table(allchunks_diamondnr_andblastntclustered_otherEukaryota_withsequence, file = paste0("taxonomy_hits_otherEukaryota_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

fasta_SAR_Eukaryotes <- read.fasta("taxonomy_hits_SAReukaryotes_list.fasta")
## rename seq.name to fullquery
fasta_SAR_Eukaryotes <- fasta_SAR_Eukaryotes %>% rename(query = seq.name)
allchunks_diamondnr_andblastntclustered_SAR_Eukaryotes_withsequence <- left_join(allchunks_diamondnr_andblastntclustered_SAR_Eukaryotes, fasta_SAR_Eukaryotes)
write.table(allchunks_diamondnr_andblastntclustered_SAR_Eukaryotes_withsequence, file = paste0("taxonomy_hits_SAReukaryotes_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# all nofishnohuman - allchunks_blastnanddiamond_hits_nofishnohuman_list.fasta
fasta_nofishnohuman <- read.fasta("taxonomy_hits_nonhost_list.fasta")
## rename seq.name to fullquery
fasta_nofishnohuman <- fasta_nofishnohuman %>% rename(query = seq.name)
allchunks_diamondnr_andblastntclustered_withsequence <- left_join(allchunks_diamondnr_andblastntclustered, fasta_nofishnohuman)
## save only 
#write.table(allchunks_diamondnr_andblastntclustered_withsequence, file = paste0("taxonomy_hits_nonhost_withsequence_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

#write.table(allchunks_diamondnr_andblastntclustered_withsequence, file = paste0("taxonomy_hits_nonhost_withsequence_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write_tsv(allchunks_diamondnr_andblastntclustered_withsequence, file = paste0("taxonomy_hits_nonhost_withsequence_mostrecent.tsv.gz"))


## then need to import into Geneious, including other fields
## the batch rename to add these to sequence name...