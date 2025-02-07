#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(scales)
library(fs)
library(taxonomizr)
library(tidyverse)
library(stringr)

parser <- ArgumentParser(description = "NR Diamond Processing")
parser$add_argument("-i", "--input")
args <- parser$parse_args()

chunk_file <- args$input

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# TODO: Handle case w/ more than 26 sub-chunks
chunk_pattern <- "chunk_\\d{3}_nonzfhum\\.\\w"
chunk_name <- stringr::str_extract(chunk_file, chunk_pattern)

single_hit <- fread(chunk_file, quote = "", header = FALSE)
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

## NOTE SECOND TIME CHANGE TO DIAMOND!!
single_hit$analysis <- "diamond"

## WILL WANT TO SEPARATELY PROCESS BLASTN THEN DIAMOND OUTPUTS, THEN AT END FULL_JOIN & FIND WHAT IS MISSING

## 60K BLASTN AND DIAMOND OUTPUTS
# (base) ➜  wikipathways head sixPRJtest_E_ntv5full
# PRJNA293022_E_NODE_1_length_9439_cov_9428.866560_g0_i0	gi|2329100573|ref|XR_007966196.1|	109293	Hippocampus zosterae	Eukaryota	75.101	494	99	5	21	3669	4150	2100	2581	7.77e-48	209
## qseqid sseqid staxid ssciname sskingdom pident length mismatch qcovs gapopen qstart qend sstart send evalue bitscore stitle
## MAY 20 - ADDING stitle - "target_title" to BOTH BLASTN & DIAMOND...TO END X17

# (base) ➜  wikipathways head sixPRJtest_E_nrdiamondfull
# PRJNA293022_E_NODE_1_length_9439_cov_9428.866560_g0_i0	XP_051788188.1	27687	Erpetoichthys calabaricus	Eukaryota	60.9	892	343	28.3	5	3450	6119	2	889	0.0	1147
## qseqid sseqid staxids sscinames sskingdoms pident length mismatch qcovhsp gapopen qstart qend sstart send evalue bitscore stitle NOTE QCOV IS SLIGHTLY DIFFERENT FOR DIAMOND!
## NOTE DIAMOND USES staxids sscinames WHICH MEANS WE HAVE TO CLEAN UP THOSE 2 COLUMNS...(I sent a github issue to Diamond author in case he can fix first)

## rename tax_division to superkingdom & tax_phylum to phylum


## RIGHT AWAY CALCULATE totalalnlen - BEFORE THE 99% MAXBITS FILTER!!
## then the 99% MAXBITS FILTER
## then CALCULATE LCA - ONLY after LCA step should we slice to get one hit per query

# we need to separate first to get length...
single_hit <- single_hit %>% separate(query, into = c("bioproject", "remove", "remove1", "node", "remove2", "length", "remove3", "coverage", "gene", "allele"), sep = "_", remove = FALSE, convert = TRUE, extra = "drop")
single_hit$remove <- NULL
single_hit$remove1 <- NULL
single_hit$remove2 <- NULL
single_hit$remove3 <- NULL


### totalalnlen calculation. JUNE 2024 - NOTE FOR NR NEED TO MULTIPLE ALNLEN BY 3 !!!
single_hit <- single_hit %>% arrange(desc(query),desc(bits))
#single_hits <- single_hit %>% group_by(query) %>% summarize(meanpident = mean(pident))
single_hitt <- single_hit %>% group_by(query) %>% summarize(sumalnlen = sum(alnlen * 3))

#single_hitu <- left_join(single_hit,single_hits)
#single_hit <- left_join(single_hitu,single_hitt)
single_hit <- left_join(single_hit,single_hitt)
## then create a new column sumperc_cov JUNE 2024 - NOTE FOR NR NEED TO MULTIPLE ALNLEN BY 3 !!!, not length
## first force as.numeric
single_hit$length <- as.numeric(single_hit$length)
single_hit$coverage <- as.numeric(single_hit$coverage)

single_hit$sumperc_cov <- single_hit$sumalnlen / single_hit$length
single_hit <- single_hit %>% relocate(sumperc_cov, .after = pident)
#####



## now 99% maxbits filter
single_hit_maxbits <- single_hit %>% group_by(query) %>% summarize(maxbits = max(bits))
single_hitss <- left_join(single_hit, single_hit_maxbits)
single_hitss$bits_percmax <- single_hitss$bits / single_hitss$maxbits
top_hits <- single_hitss %>% dplyr::filter(bits_percmax > 0.989999)

top_hits$rel_abundance <- (top_hits$length * top_hits$coverage * 0.01)
top_hits <- top_hits %>% relocate(rel_abundance, .after = coverage)


# rm(single_hit)
# rm(single_hitss)
# rm(single_hit_maxbits)

## some cleaning before LCA

## also removing artificial sequences, other, etc. at step A now
top_hits <- top_hits %>% filter(!grepl("N/A", top_hits$taxname))
top_hits <- top_hits %>% filter(!grepl("vector", top_hits$taxname))
top_hits <- top_hits %>% filter(!grepl("Vector", top_hits$taxname))
top_hits <- top_hits %>% filter(!grepl("synthetic", top_hits$taxname))
#top_hits <- top_hits %>% filter(!grepl(";", top_hits$taxid))

## INSTEAD OF REMOVING ALL ; , REMOVE AFTER FIRST INSTANCE OF ; IN TAXID & AFTER 0 IN TAXNAME
# top_hitscheck <- top_hits %>% filter(grepl("0", top_hits$taxname))
# top_hits %>% filter(grepl(";", top_hits$taxid))
# top_hits0 <- top_hits
top_hits$taxid <- gsub(";.*", "", top_hits$taxid)
top_hits$taxname <- gsub("0.*", "", top_hits$taxname)

##need to convert taxid to as.numeric
top_hits$taxid <- as.numeric(top_hits$taxid)

## often way too many taxids, so extract only unique
n_distinct(top_hits$taxid)
taxid_unique <- unique(top_hits$taxid)

## also clade for SAR???
taxa <- getTaxonomy(taxid_unique, '/tax.sqlite', desiredTaxa = c("superkingdom", "clade", "kingdom", "phylum", "class", "order", "family", "genus","species"))
taxa.df <- as.data.frame(taxa)
taxa.df <- as_tibble(taxa, rownames = "taxid")
taxa.df$taxid <- as.numeric(taxa.df$taxid)

## now need to unlabel rows and name column
top_hitswithtaxa <- left_join(top_hits,taxa.df)

rm(top_hits)
top_hits <- top_hitswithtaxa

########### END TAXONOMIZR CODE ###########

top_hits <- top_hits %>% select(-any_of(c("skingdom", "tax_kingdom")))
top_hits <- top_hits %>% rename("tax_superkingdom" = "superkingdom") %>% rename("tax_clade" = "clade") %>% rename("tax_kingdom" = "kingdom") %>% rename("tax_phylum" = "phylum") %>%rename("tax_class" = "class") %>% rename("tax_order" = "order") %>% rename("tax_family" = "family") %>% rename("tax_genus" = "genus") %>% rename("tax_species" = "species")


## now using mmseqs2 code

## adding curated last taxoncategory column ;-_Sar
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
top_hits <- top_hits %>% mutate(across(taxoncategory, ~ replace_na(., "other_Eukaryota")))

## note there are both zeroes & NAs in this...
# top_hits$tax_kingdom <- gsub("0", "", top_hits$tax_kingdom)
# top_hits$phylum <- gsub("0", "", top_hits$phylum)
## weird one-ff if tax_kingdom is all zeros
top_hits$tax_kingdom <- as.character(top_hits$tax_kingdom)
top_hits <- top_hits %>%
mutate(tax_kingdom = na_if(tax_kingdom, "0"))
top_hits <- top_hits %>%
mutate(tax_phylum = na_if(tax_phylum, "0"))

## HERE USE BLASTN FIRST, THEN DIAMOND NEXT
## ALSO FOR DIAMOND, NEED TO APPEND "_diamond" to all columns, then do a full join
## ACTUALLY APPEND WITH _NR - ALSO APPEND WITH _NT FOR BLASTN
# ## Adding suffixes
# df %>% rename_with( ~ paste0(.x, "_diamond"))
# df %>% rename_with( ~ paste0(.x, "_NT"))
# df %>% rename_with( ~ paste0(.x, "_NR"))

thresholded_hit_nofishnomammals0 <- top_hits %>% filter(!grepl("Actinopteri", top_hits$tax_class))
thresholded_hit_nofishnomammals0 <- thresholded_hit_nofishnomammals0 %>% filter(!grepl("Primates", thresholded_hit_nofishnomammals0$tax_order))

# write.table(thresholded_hit_nofishnomammals_diamond, file = paste0(currentrunname,"diamond_hits_test2",".tab"), sep = "\t", row.names = FALSE, quote = FALSE)

#################### CODE TO PULL LCA FROM TAXONOMY ##############################

# MAY 2024 - NEED TO UPDATE THIS, GROUP BY QUERY,TAXID - mutate (new column) that is n_distinct (query,taxid)
#mutate(mass_norm = mass / mean(mass, na.rm = TRUE))
## also copy this to first pass!

## first arrange thresholded_hit_nofishnomammals0
thresholded_hit_nofishnomammals0 <- thresholded_hit_nofishnomammals0 %>% arrange(bioproject,desc(length),desc(node),desc(bits))
## now slice by query and taxid
thresholded_hit_nofishnomammals7 <- thresholded_hit_nofishnomammals0 %>% group_by(query,taxid) %>% slice_head(n = 1) %>% ungroup()
# # blast_homodanio_PRJNA381153 <- blast_homodanio_PRJNA381153 %>% group_by(query,target) %>% slice_head(n = 1) %>% ungroup()

thresholded_hit_nofish_dups <- thresholded_hit_nofishnomammals7 %>% group_by(query) %>% add_tally() %>% dplyr::filter(n > 1) %>% ungroup()
thresholded_hit_nofish_rest <- thresholded_hit_nofishnomammals7 %>% group_by(query) %>% add_tally() %>% dplyr::filter(n == 1) %>% ungroup()
# 

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




### 5/17 - add more from below (there is an if/else, can just run the main set of code)

## lots of steps here...NOTE THE na.rm = true helps in some cases, it is much worse in others, skipping for now
#   examining_dups1 <- thresholded_hit_nofish_dups %>% group_by(query) %>% summarise(across(superkingdom, ~n_distinct(., na.rm = TRUE)))
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

## works once, but we want nested https://www.statology.org/nested-ifelse-in-r/
#thresholded_hit_nofish_dups2$taxname_new <- ifelse((thresholded_hit_nofish_dups2$taxlineage2_n > 1), thresholded_hit_nofish_dups2$taxlineage1, 0)


thresholded_hit_nofish_dups2$taxname_lca <- ifelse((thresholded_hit_nofish_dups2$taxlineage1_n > 1), "root", 
                                                 ifelse((thresholded_hit_nofish_dups2$taxlineage2_n > 1), thresholded_hit_nofish_dups2$tax_superkingdom, 
                                                        ifelse((thresholded_hit_nofish_dups2$taxlineage3_n > 1), thresholded_hit_nofish_dups2$tax_clade,
                                                               ifelse((thresholded_hit_nofish_dups2$taxlineage4_n > 1), thresholded_hit_nofish_dups2$tax_kingdom,
                                                                      ifelse((thresholded_hit_nofish_dups2$taxlineage5_n > 1), thresholded_hit_nofish_dups2$tax_phylum,
                                                                             ifelse((thresholded_hit_nofish_dups2$taxlineage6_n > 1), thresholded_hit_nofish_dups2$tax_class,
                                                                                    ifelse((thresholded_hit_nofish_dups2$taxlineage7_n > 1), thresholded_hit_nofish_dups2$tax_order,
                                                                                           ifelse((thresholded_hit_nofish_dups2$taxlineage8_n > 1), thresholded_hit_nofish_dups2$tax_family, thresholded_hit_nofish_dups2$tax_genus))))))))

thresholded_hit_nofish_dups2 <- thresholded_hit_nofish_dups2 %>% relocate(taxname_lca, .after = taxname)


## not quite right, because if taxon is missing may add a zero to the taxname_lca
#thresholded_hit_nofish_dups2 <- thresholded_hit_nofish_dups2 %>% mutate(across(taxname_lca, ~ replace_na(., superkingdom)))

thresholded_hit_nofish_dups2$taxname_lca <- ifelse((is.na(thresholded_hit_nofish_dups2$taxname_lca) == TRUE), thresholded_hit_nofish_dups2$tax_superkingdom, thresholded_hit_nofish_dups2$taxname_lca)

## checking to see what reads are completely unresolved to root...
#thresholded_hit_nofish_rootcheck <- thresholded_hit_nofish_dups2 %>% filter(taxname_lca == "root")
#unique(thresholded_hit_nofish_dups2$taxname_lca)

## extract just query and taxname new, then rename that...then add to thresholded_hit_nofish_rest to make new thresholded_hit_nofish
thresholded_hit_nofish_dups3 <- thresholded_hit_nofish_dups2 %>% select(query,taxname_lca)
#  rm(thresholded_hit_nofish_dups2)
## changed from thresholded_hit_nofishnomammals0 to thresholded_hit_nofishnomammals7
thresholded_hit_nofish_dups4 <- thresholded_hit_nofishnomammals7 %>% group_by(query) %>% add_tally() %>% dplyr::filter(n > 1)
thresholded_hit_nofish_dups3 <- thresholded_hit_nofish_dups3 %>% group_by(query) %>% slice(1)
thresholded_hit_nofish_dups5 <- left_join(thresholded_hit_nofish_dups4, thresholded_hit_nofish_dups3) %>% 
select(-taxname) %>% group_by(query) %>% slice(1) %>% relocate(taxname_lca, .after = taxid) %>% select(-taxoncategory)
rm(thresholded_hit_nofish_dups3)
rm(thresholded_hit_nofish_dups4)


## also for dups5, want to refresh the taxon category
## but we want to keep taxname_lca - maybe keep taxoncategory2
#thresholded_hit_nofish_dups5$taxoncategory <- thresholded_hit_nofish_dups5$taxname_lca
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


thresholded_hit_nofishnomammals <- thresholded_hit_nofishnomammals %>% mutate(across(taxoncategory2, ~ replace_na(., "other_Eukaryota")))
rm(thresholded_hit_nofish_dups5)
rm(thresholded_hit_nofish_rest)
rm(examining_dups_all)

## UPDATE 9/7 - REMOVE TARGET_TITLE & ENSURE ORDER OF 4 COLUMNS
thresholded_hit_nofishnomammals <- thresholded_hit_nofishnomammals %>% relocate(target, .after = rel_abundance) %>% relocate(bits, .after = target) %>% relocate(taxid, .after = bits) %>% relocate(taxname_lca, .after = taxid)
if("target_title" %in% colnames(thresholded_hit_nofishnomammals) == FALSE) {
thresholded_hit_nofishnomammals$target_title <- "N/A"
thresholded_hit_nofishnomammals <- thresholded_hit_nofishnomammals %>% relocate(target_title, .before = taxoncategory2)
}

## rename columns 8-36 by appending _NR
thresholded_hit_nofishnomammals_diamond <- thresholded_hit_nofishnomammals %>% rename_with(~paste0(.x, "_NR"), c(target:taxoncategory2))


## INSTEAD USE CURRENT CHUNKNAME
## chunk names
currentchunkname <- chunk_name

## SAVE THE OUTPUTS TO CSV, THEN STREAMLINE CODE FOR R ONLY...
write.table(thresholded_hit_nofishnomammals_diamond, file = paste0(currentchunkname, ".diamond_hits.tab"), sep = "\t", row.names = FALSE, quote = FALSE)
