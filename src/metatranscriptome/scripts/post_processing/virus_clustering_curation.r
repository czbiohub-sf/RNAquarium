library(data.table)
library(scales)
library(fs)
library(taxonomizr)
library(tidyverse)

library(phylotools)
library(treemapify)

###
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
## Mode2 ignores NA
Mode2 <- function(x) {
  ux <- na.omit(unique(x) )
  ux[which.max(tabulate(match(x, ux)))]
}

workingpath <- getwd()
workingpathdash <- str_c(workingpath, "/")
diamondpath <- str_c(workingpath, "/nr_diamond")
blastpath <- str_c(workingpath, "/nt_blast")

## saving all outputs to RNAquarium_outputs/virus_outputs
outpath <- str_c(workingpath, "/RNAquarium_outputs")
#setwd(outpath)

outpathvirus <- str_c(outpath, "/virus_outputs")
setwd(outpathvirus)


## for clustering, we want to create new cluster directories for final steps - save for next step
#clusters_dir_name <- paste0("clusters_", format(Sys.Date(), "%m%d%y"))
#dir.create(clusters_dir_name)

### fasta-writing functions
writetoFasta <- function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"read"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq.text"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}
## optimized
writetoFastafaster <- function(data, filename) {
  # Construct the FASTA lines in a vectorized way
  fastaLines <- paste0(">", data[["read"]], "\n", data[["seq.text"]])
  
  # Write directly to the file
  writeLines(fastaLines, filename)
}


####################################
## after mmseqs2 easy-cluster R commands to curate the results


#library(phylotools) taxonomy_hits_viruses_mostrecent_withtargets.fasta taxonomy_hits_viruses_easycluster1_cluster.tsv
clusters_viruses_fasta <- read.fasta("taxonomy_hits_viruses_mostrecent_withtargets.fasta")
clusters_viruses_fasta <- clusters_viruses_fasta %>% rename(read = seq.name)

### for targets, need to remove spaces - now done in previous step
#clusters_viruses_fasta$read <- gsub(' ','_',clusters_viruses_fasta$read)

clusters_viruses <- read_tsv("taxonomy_hits_viruses_easycluster1_cluster.tsv", col_names = FALSE)
clusters_viruses <- clusters_viruses %>% rename(cluster = X1) %>% rename(read = X2)

## then group by cluster
clusters_viruses_tally0 <- clusters_viruses %>% group_by(cluster) %>% add_tally()




## first batch of largest (usually none here, but check)
# clusters_viruses_tally <- clusters_viruses_tally0 %>% dplyr::filter(n > 5000)
# ## this is as single larger cluster of caudoviricetes, 9k!

## part 1
#clusters_viruses_tally1 <- clusters_viruses_tally0 %>% dplyr::filter(n > 2500)
clusters_viruses_tally <- clusters_viruses_tally0 %>% dplyr::filter(n > 9)

n_distinct(clusters_viruses_tally$cluster)

## part 2 all between 1 & 9!!
clusters_viruses_tally2 <- clusters_viruses_tally0 %>% dplyr::filter(n < 10)
clusters_viruses_tally2 <- clusters_viruses_tally2 %>% dplyr::filter(n > 1)
## also singletons
clusters_viruses_singletons <- clusters_viruses_tally0 %>% dplyr::filter(n == 1)


### some summaries then save for combining with allchunks_blastnanddiamond_hits_viruses_withsequence_mostrecent.tsv
n_distinct(clusters_viruses_tally$cluster)
n_distinct(clusters_viruses_tally2$cluster)
n_distinct(clusters_viruses_singletons$cluster)


clusters_viruses_summary <- inner_join(clusters_viruses_fasta,clusters_viruses_tally0)
## also can arrange by length??
clusters_viruses_summary$length <- nchar(clusters_viruses_summary$seq.text)
clusters_viruses_summary <- clusters_viruses_summary %>% select(read,cluster,n,length)

clusters_viruses_summary <- clusters_viruses_summary %>% arrange(desc(n),cluster,desc(length)) %>% group_by(cluster)

clusterorder <- clusters_viruses_summary %>% arrange(desc(n),cluster,desc(length))
clusterorder2 <- unique(clusterorder$cluster)
#clusterorder <- clusterorder[["cluster"]]
clusters_viruses_summary$cluster <- factor(clusters_viruses_summary$cluster, levels = clusterorder2, ordered = TRUE)

clusters_viruses_summary <- clusters_viruses_summary %>% slice_head(n = 1)
#clusters_viruses_summary <- clusters_viruses_summary %>% dplyr::filter(n == 1 & grepl('_NODE', read)) dplyr::filter(!grepl('phage|Phage', taxname_lca))
clusters_viruses_summary <- clusters_viruses_summary %>% dplyr::filter((grepl('_NODE', read) & n == 1) | n > 1)
clusters_viruses_summary
write.table(clusters_viruses_summary, file = paste0("taxonomy_hits_viruses_clustersummary",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)              
write.table(clusters_viruses_summary, file = paste0("taxonomy_hits_viruses_clustersummary_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)              


## now actually combine taxonomy_hits_viruses_mostrecent_withtargets.fasta

# ### now add "cluster_size_n" to start of each cluster in this file, then add to thresholded_hit_viruses_NTNR_withsequence
# clusters_viruses_tally0$cluster0 <- "cluster_size"
# #clusters_viruses_tally0$cluster1 <- "size"
# 
# ## unite
# clusters_viruses_tally0 <- clusters_viruses_tally0 %>% unite(clusterfull, cluster0, n, cluster, sep = "_", remove = FALSE)
# clusters_viruses_tally0$cluster0 <- NULL
# clusters_viruses_tally0$cluster <- NULL
# 
# clusters_viruses_tally0 <- clusters_viruses_tally0 %>% rename(query = read)
# clusters_viruses_tally0 <- clusters_viruses_tally0 %>% rename(clustersize = n)
# 
# ## now left join
#allchunks_diamondnr_andblastntclustered_viruses_withsequence0 <- allchunks_diamondnr_andblastntclustered_viruses_withsequence
allchunks_diamondnr_andblastntclustered_viruses_withsequence <- read_tsv("taxonomy_hits_viruses_withsequence_mostrecent.tsv")

targetNT_fasta <- clusters_viruses_fasta %>% dplyr::filter(!grepl('P_NODE_', read))
targetNT_fasta <- targetNT_fasta %>% dplyr::filter(!grepl('S_NODE_', read))

## also create some columns then join to allchunks_diamondnr_andblastntclustered_viruses_withsequence
targetNT_fasta <- targetNT_fasta %>% dplyr::rename(query = read)
targetNT_fasta$bioproject <- "NTtarget"
targetNT_fasta$length <- nchar(targetNT_fasta$seq.text)
targetNT_fasta$analysis_used <- "NTtarget"
targetNT_fasta$taxname_lca_NTorNR <- "NTtarget"
targetNT_fasta$taxoncategory_NTorNR <- "NTtarget"
targetNT_fasta$taxoncategorysimple_NTorNR <- "NTtarget"
#targetNT_fasta$bits_NTorNR <- "NTtarget"
#targetNT_fasta$evalue_NTorNR <- "NTtarget"
targetNT_fasta$viruscategory_NTorNR <- "NTtarget"
targetNT_fasta$viruscategorysimple_NTorNR <- "NTtarget"
targetNT_fasta <- targetNT_fasta %>% relocate(seq.text, .after = viruscategorysimple_NTorNR)

#allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence0 <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence
## aug 2025 - should be bind_rows not full_join
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence <- bind_rows(allchunks_diamondnr_andblastntclustered_viruses_withsequence,targetNT_fasta)
#allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence <- full_join(allchunks_diamondnr_andblastntclustered_viruses_withsequence,targetNT_fasta)


## now we need to make clusters_viruses_tally0 read look like query in allchunks_diamondnr_andblastntclustered_viruses_withsequence
#PRJ550047_P_NODE_615734_length_234_cov_17.354037_g612699_i0|Caudoviricetes_sp.|NTclustered|2.95e-115|
# to
#PRJDB4470_P_NODE_67321_length_254_cov_2.429268_g63478_i0
## need a separate()

clusters_viruses_tally0_tomerge <- clusters_viruses_tally0 %>% separate_wider_delim(read, delim = "|", names = c("query"), too_many = "drop")


allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2 <- left_join(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence,clusters_viruses_tally0_tomerge)


allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2 <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2 %>% relocate(cluster, .after = query)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2 <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2 %>% relocate(n, .after = cluster)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2 <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2 %>% rename(clustersize = n)
## export to merge & look...
#write.table(clusters_viruses_tally0_tomerge, file = paste0("clusters_viruses_tally0_tomerge_testing",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## missing - look identical???
# PRJNA572587_P_NODE_66432_length_77_cov_11.507937_g65451_i0
# PRJNA572587_P_NODE_66432_length_77_cov_11.507937_g65451_i0
# PRJ572587_P_NODE_66432_length_77_cov_11.507937_g65451_i0
## note they are mssing some of the start!! IT IS BECAUSE OF A |NA REPLACE IN VIRUS_CURATION STEP



## also need to pull in target IDs before left_joining with allchunks_diamondnr_andblastntclustered_viruses_withsequence
## again always have a generic version saved for after clustering
#write.table(allchunks_diamondnr_andblastntclustered_viruses_withsequence, file = paste0("allchunks_blastnanddiamond_hits_viruses_withsequence_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(thresholded_hit_viruses_NTNR_withsequenceandclusters, file = paste0("RNaquarium_allchunks_blastnanddiamond_hits_viruses_nonphagelist_withsequenceandclusters",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


## remove singleton NTtargets
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2 <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2 %>% dplyr::filter((grepl('_NODE', query) & clustersize == 1) | clustersize > 1)


## create bioprojectsbycluster
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2a <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2 %>% group_by(cluster) %>% summarize(bioprojectsbycluster=n_distinct(bioproject))
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- left_join(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2,allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2a)
# rm(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2a)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2 %>% group_by(cluster) %>% mutate(bioprojectsbycluster=n_distinct(bioproject))
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(bioprojectsbycluster, .before = node) %>% ungroup()

## different filter - this also will remove clusters of two that are only nttargets
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% dplyr::filter(!(bioproject == "NTtarget" & bioprojectsbycluster == 1))

## after filter - recalculate bioprojectsbycluster to NOT include NTtarget as bioproject!
## instead of filtering, convert to NA and use n_distinct(..., na.rm = TRUE)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% mutate(bioproject = na_if(bioproject, "NTtarget"))
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% group_by(cluster) %>% mutate(bioprojectsbycluster2 = n_distinct(bioproject, na.rm = TRUE)) %>% ungroup()
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(bioprojectsbycluster2, .after = bioprojectsbycluster) %>% ungroup()
## now add back NTtarget using replace_na() replace_na(list(cpart2 = ""))
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% replace_na(list(bioproject = "NTtarget"))
## also remove bioprojectsbycluster & rename bioprojectsbycluster2
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% select(-bioprojectsbycluster) %>% rename(bioprojectsbycluster = bioprojectsbycluster2)


### now to make a full name (for virus counts, etc) also move columns. "cluster"
# ### now add "cluster_size_n" to start of each cluster in this file, then add to thresholded_hit_viruses_NTNR_withsequence
# clusters_viruses_tally0$cluster0 <- "cluster_size"
# #clusters_viruses_tally0$cluster1 <- "size"

## FEB 2025 - SHORTENING QUERY AND USING CREATING A SHORTER FULLQUERY FIELD AS WELL
## first shorten query by splitting back up, only keeping prj and node only separate() then unite()
## already have bioproject & node, can do this with just a unite
#PRJxxx_CONTIG_(node)|
  

## need an if/then for the NTtargets
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$drop <- "CONTIG"

## create NTtargets if else field
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% mutate(is_target = str_detect(bioproject, "NTtarget")) 

## if NTtargets, 
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% mutate(query2 = substr(query, 1, 100))
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% unite(shortqueryNTyes, bioproject, query2, sep = "_", remove = FALSE)
## if not NTtargets:
#allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$drop <- "CONTIG"
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% unite(shortqueryNTno, bioproject, drop, node, sep = "_", remove = FALSE)
## if else here
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% mutate(shortquery = if_else(is_target,shortqueryNTyes,shortqueryNTno)) 

allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$drop <- NULL
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$query2 <- NULL
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$is_target <- NULL
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$shortqueryNTyes <- NULL
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$shortqueryNTno <- NULL

allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(shortquery, .after = query)

allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$shortquery <- gsub('_MAG_TPA_asm_','_',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$shortquery)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$shortquery <- gsub('_TPA_asm_','_',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$shortquery)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$shortquery <- gsub('_MAG_','_',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$shortquery)

## similar if else for LCA - only run NO JUST MERGE WITH LCA, IT WILL BE "NTCLUSTER" taxname_lca_NTorNR
## then combine to make full query as before but fewer terms
##  combine with lca

allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% unite(fullquery, shortquery, taxname_lca_NTorNR, sep = "|", remove = FALSE)

## add in here a check for rare cases where shortquery name is duplicated (when NODE_P & NODE_S are otherwise identical)...

## create a new column: nshortquery - group_by(shortquery) %>% mutate(nshortquery=n()) %>% relocate(nshortquery, .after = shortquery)
## also create new columns by separate() query by '_' - use second group and can drop rest
## separate_wider_delim(query, delim = "_", names = c("drop", "query_sorp"), too_many = "drop", too_few = "align_start", cols_remove = FALSE)
## for each row: if nshortquery > 1, then take shortquery column and unite() with query_sorp %>% unite(shortquery2, shortquery, query_sorp, sep = "_")
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% group_by(shortquery) %>% mutate(nshortquery=n()) %>% relocate(nshortquery, .after = shortquery) %>% ungroup()


allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b_dups <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% dplyr::filter(nshortquery > 1)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b_nodups <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% dplyr::filter(nshortquery == 1)

allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b_dups <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b_dups %>% separate_wider_delim(query, delim = "_", names = c("drop", "query_sorp"), too_many = "drop", too_few = "align_start", cols_remove = FALSE)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b_dups <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b_dups %>% unite(shortquery2, shortquery, query_sorp, sep = "_")

## rename shortquery2 to shortquery, also grep and change _P and _S to p & s - actually either gsub or str_replace_all mutate(taxname_virus = str_replace_all(taxname_virus, c("Sinsheimervirus" = "phiX phage")))
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b_dups <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b_dups %>% select(-drop) %>% rename(shortquery = shortquery2)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b_dups <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b_dups %>% mutate(shortquery = str_replace_all(shortquery, c("_P" = "p"))) %>%  mutate(shortquery = str_replace_all(shortquery, c("_S" = "s")))

## join again or bind_rows, finally drop nshortquery
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- bind_rows(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b_nodups,allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b_dups)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% select(-nshortquery) %>% ungroup()
rm(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b_dups)
rm(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b_nodups)


## next re-name clustername - again shortening, FIND THE MOST COMMON LCA IN THE CLUSTER
## check to see if cluster LCA + n + bioprojects will be unique!!  
# ## Mode2 ignores NA
# Mode2 <- function(x) {
#   ux <- na.omit(unique(x) )
#   ux[which.max(tabulate(match(x, ux)))]
# }

## group_by() then mode()
## need to first duplicate taxname_lca_NTorNR but change NTtarget to NA
## this isn't working well, we instead want to create another version of taxname_lca_NTorNR2
## in this case if is_NTtarget, then make taxname_lca_NTorNR2 = shortquery
## otherwise taxname_lca_NTorNR2 = taxname_lca_NTorNR
## then we don't need to worry about NAs and can use Node again

## create NTtargets if else field
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% mutate(is_target = str_detect(bioproject, "NTtarget")) 

allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% mutate(taxname_lca_NTorNR2 = if_else(is_target,shortquery,taxname_lca_NTorNR)) 

#allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$taxname_lca_NTorNR2 <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$taxname_lca_NTorNR
#allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% mutate(taxname_lca_NTorNR2 = na_if(taxname_lca_NTorNR2, "NTtarget"))

#allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b_clustermode <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% group_by(cluster) %>% mutate(clusterLCA = Mode2(taxname_lca_NTorNR))
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% group_by(cluster) %>% mutate(clusterLCA = Mode(taxname_lca_NTorNR2))
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(clusterLCA, .after = cluster) %>% ungroup()

# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(clusterLCA, .after = cluster)
# ## then a join
# # WE ALSO WANT TO ADD A SIMPLE ID# TO EACH CLUSTER BY ITS CLUSTERLCA
#  
# ## first arrange by desc
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% arrange(clusterLCA,desc(clustersize),node)
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c %>% group_by(clusterLCA) %>% mutate(clusterLCAid = as.numeric(as.factor(cluster))) %>% relocate(clusterLCAid, .after = clusterLCA) %>% ungroup()
# # 
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% group_by(clusterLCA) %>% arrange(clusterLCA, desc(clustersize), node, .by_group = TRUE) %>% mutate(clusterLCAid = as.numeric(as.factor(cluster))) %>% relocate(clusterLCAid, .after = clusterLCA) %>% ungroup()
# # 
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% arrange(clusterLCA,desc(clustersize),node)
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% group_by(clusterLCA) %>% arrange(desc(clustersize),node, .by_group = TRUE) %>% mutate(clusterLCAid = as.numeric(as.factor(cluster))) %>% relocate(clusterLCAid, .after = clusterLCA) %>% ungroup()
# # 
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% arrange(clusterLCA,desc(clustersize),node)
# 
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% group_by(clusterLCA) %>% arrange(clusterLCA, desc(clustersize), node, .by_group = TRUE) %>% mutate(clusterLCAid = dense_rank(cluster)) %>% relocate(clusterLCAid, .after = clusterLCA) %>% ungroup()
# # 
# # #allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% group_by(clusterLCA) %>% arrange(clusterLCA, desc(clustersize), node, .by_group = TRUE) %>% mutate(clusterLCAid = as.numeric(as.factor(cluster))) %>% relocate(clusterLCAid, .after = clusterLCA) %>% ungroup()
# # #allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% arrange(desc(clustersize)) %>% group_by(clusterLCA) %>% mutate(clusterLCAid = as.numeric(as.factor(cluster))) %>% relocate(clusterLCAid, .after = clusterLCA) %>% ungroup()
# # 
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% group_by(clusterLCA) %>% arrange(node, desc(clustersize), .by_group = TRUE) %>% mutate(clusterLCAid = dense_rank(cluster)) %>% relocate(clusterLCAid, .after = clusterLCA) %>% ungroup()
# 
# rm(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c)
# #allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% arrange(desc(clustersize),clusterLCA,node)
# #allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c %>% group_by(clusterLCA) %>% mutate(clusterLCAid = as.numeric(as.factor(cluster))) %>% relocate(clusterLCAid, .after = clusterLCA) %>% ungroup()
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c %>% group_by(clusterLCA) %>% mutate(clusterLCAid = dense_rank(cluster)) %>% relocate(clusterLCAid, .after = clusterLCA) %>% ungroup()
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c %>% group_by(clusterLCA) %>% arrange(desc(clustersize),node) %>% mutate(clusterLCAid = dense_rank(cluster)) %>% relocate(clusterLCAid, .after = clusterLCA) %>% ungroup()
# 
# 
# 
# ## make new df with factor for correct order?
# # idorder <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% arrange(clusterLCA, desc(clustersize), node)
# # idorder2 <- unique(idorder$cluster)
# # #clusterorder <- clusterorder[["cluster"]]
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c$cluster <- factor(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$cluster, levels = idorder2, ordered = TRUE)
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c %>% mutate(clusterLCAid = dense_rank(cluster)) %>% relocate(clusterLCAid, .after = clusterLCA) %>% ungroup()
# 
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b
# # idorder <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c %>% arrange(clusterLCA, desc(clustersize), node)
# # idorder2 <- unique(idorder$cluster)
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c$cluster <- factor(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c$cluster, levels = idorder2, ordered = TRUE)
# # #allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c %>% mutate(clusterLCAid = dense_rank(cluster)) %>% relocate(clusterLCAid, .after = clusterLCA) %>% ungroup()
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2c %>% mutate(clusterLCAid = as.numeric(as.factor(cluster))) %>% relocate(clusterLCAid, .after = clusterLCA) %>% ungroup()
# 
# ## finally combine fullquery and clusterLCA also clustersize & bioprojectsbycluster to make fullquery_clusterLCA
# # PRJxxx_CONTIG_(node)|LCA|CLUSTER_(clusterLCA)_n|(numbioprj)
# ## don't need code to distiguish clustername anymore...    
# 
# # ## first create clustername depending on if cluster is a contig or a target
# # ## fix - just make clustername1 & clustername2 and pick one or the other
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% mutate(clustername2 = substr(cluster, 1, 100))
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% separate_wider_delim(cluster, delim = "|", names = c("cpart1", "cpart2"), too_many = "drop", too_few = "align_start", cols_remove = FALSE)
# # ## reverse order, have lca name first!!
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% unite(clustername1, cpart2, cpart1, sep = "|")
# # 
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% mutate(has_pipe = str_detect(cluster, "\\|")) 
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% mutate(clustername = if_else(has_pipe,clustername1,clustername2)) 
# # 
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% select(-has_pipe) %>% select(-clustername1) %>% select(-clustername2)
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(clustername, .before = clustersize)
# # ## remove _MAG_TPA_asm_ etc clusters_viruses_fasta_smaller$read <- gsub('NA$','',clusters_viruses_fasta_smaller$read)
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$clustername <- gsub('_MAG_TPA_asm_','_',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$clustername)
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$clustername <- gsub('_TPA_asm_','_',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$clustername)
# # allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$clustername <- gsub('_MAG_','_',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$clustername)
# 
# ## now start to combine - first by re-making larger fullquery name fullquery and clusterLCA also clustersize & bioprojectsbycluster to make fullquery_clusterLCA
# #allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% unite(fullquery, query, taxname_lca_NTorNR, analysis_used, evalue_NTorNR, lowcoverage_flag, sep = "|", remove = FALSE)
# 
# 

# group_by undoes the arrange, no obvious way to do this using tidy, so using base R per https://stackoverflow.com/questions/64024505/how-do-i-restart-cur-group-id-in-r
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% arrange(desc(clustersize),node)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$clusterLCAid <- with(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b, ave(paste(cluster), clusterLCA, FUN = function(x) match(x, unique(x))))
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(clusterLCAid, .after = clusterLCA) %>% ungroup()

## HERE UPDATE TO ALSO ADD clusterLCAid MUCH LIKE THE CONTIG NAMING...

allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$drop <- "CLUSTER"
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% unite(fullquery_clusterLCA, fullquery, drop, clusterLCA, clustersize, bioprojectsbycluster, sep = "|", remove = FALSE)
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$fullquery_clusterLCA <- gsub('CLUSTER\\|','CLUSTER_',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$fullquery_clusterLCA)

#allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% unite(fullquery_clusterLCA, fullquery, drop, clusterLCAid, clusterLCA, clustersize, bioprojectsbycluster, sep = "|", remove = FALSE)
## need to do in three batches, last not using |, so we can replace...
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% unite(fullquery_clusterLCAparta, fullquery, drop, clusterLCAid, sep = "|", remove = FALSE)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% unite(fullquery_clusterLCApartb, clusterLCA, clustersize, bioprojectsbycluster, sep = "|", remove = FALSE)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% unite(fullquery_clusterLCA, fullquery_clusterLCAparta, fullquery_clusterLCApartb, sep = "?", remove = TRUE)

allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$fullquery_clusterLCA <- gsub('CLUSTER\\|','CLUSTER',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$fullquery_clusterLCA)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$fullquery_clusterLCA <- gsub('\\?','_',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$fullquery_clusterLCA)

# PRJxxx_CONTIG_(node)|LCA|CLUSTER_(clusterLCA)_n|(numbioprj)


## CHECK MAY NOT NEED THESE
# ## change |NTtarget|NTtarget|NA|NA to |NTtarget 'NTtarget\\|NTtarget\\|NA\\|NA$','NTtarget'
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$fullquery <- gsub('NTtarget\\|NTtarget\\|NA\\|NA$','NTtarget',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$fullquery)
# ## also remove |NA
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$fullquery <- gsub('\\|NA$','',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$fullquery)
# ## and spaces! 
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$fullquery <- gsub(' ','_',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$fullquery)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$fullquery_clusterLCA <- gsub(' ','_',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$fullquery_clusterLCA)

## move columns
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(query, .before = seq.text)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(cluster, .before = seq.text)


allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(tax_superkingdom_NTorNR, .before = taxname_lca_NTclustered)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(tax_clade_NTorNR, .before = taxname_lca_NTclustered)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(tax_kingdom_NTorNR, .before = taxname_lca_NTclustered)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(tax_phylum_NTorNR, .before = taxname_lca_NTclustered)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(tax_class_NTorNR, .before = taxname_lca_NTclustered)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(tax_order_NTorNR, .before = taxname_lca_NTclustered)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(tax_family_NTorNR, .before = taxname_lca_NTclustered)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(tax_genus_NTorNR, .before = taxname_lca_NTclustered)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(tax_species_NTorNR, .before = taxname_lca_NTclustered)

allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$is_target <- NULL
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$drop <- NULL
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$taxname_lca_NTorNR2 <- NULL

allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(cluster, .after = clusterLCA)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(query, .after = analysis_used)

## also get new counts and bioprojects for clusterLCA!!  #allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(clustername) %>% mutate(newcount=n())
## clusterLCA_contigcount & clusterLCA_bioprojectcount
# ## note save these counts for 'interesting' set
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% group_by(clusterLCA) %>% mutate(clusterLCA_bioprojectcount=n_distinct(bioproject))
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% group_by(clusterLCA) %>% add_tally() %>% rename(clusterLCA_contigcount = n)
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(clusterLCA_bioprojectcount, .after = clusterLCA)
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(clusterLCA_contigcount, .after = clusterLCA)


## more relocates? - note we may want to remove some of these before saving? or save two versions?
#allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$shortquery <- NULL
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$analysis_NTorNR <- NULL
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(bioproject, .before = node)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$clusterLCAid <- NULL


### making a simpler cluster name - shortclustername, very similar to shortcontigname
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% separate_wider_delim(fullquery_clusterLCA, delim = "|", names = c("drop1", "drop2", "conpart1", "conpart2", "conpart3"), too_many = "drop", too_few = "align_start", cols_remove = FALSE)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% unite(shortclustername, conpart1, conpart2, conpart3, sep = "|", remove = TRUE)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% relocate(shortclustername, .before = cluster)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$drop1 <- NULL
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b$drop2 <- NULL


## moving save to after subsetting interesting - want to rename all columns after code is complete
# ## save
# write.table(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b, file = paste0("allchunks_blastnanddiamond_hits_viruses_withsequenceandclusters_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
# ## also always have a generic version saved NOTE THIS WILL CHANGE, BUT IS USED FOR PART 3 AND ALSO VIRUS CURATION...
# write.table(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b, file = paste0("allchunks_blastnanddiamond_hits_viruses_withsequenceandclusters_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


## GET INTERESTING
## at end make an interesting subset ONLY NOW - at this point remove any row with phage Caudo Microviridae Carnation_latent Porcine_bastro 1_XXX in either read or cluster column!!
## also for 'interesting'
## remove phage Caudoviricetes Microviridae Carnation_latent Porcine_bastrovirus
## also 1_XXX

## REMOVE FROM BOTH CLUSTER & fullquery_clusterLCA
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% dplyr::filter(viruscategorysimple_NTorNR != "Unresolved Viruses")
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(viruscategorysimple_NTorNR != "Adapter")
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(viruscategorysimple_NTorNR != "Phage")
## more stringent filter here grepl Carnation_latent_virus filter(!(
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Caudoviricetes', clustername))
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('phage|Phage', clustername))
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Carnation_latent', clustername))
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Microviridae', clustername))
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Porcine_bastrovirus', clustername))
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('1_XXX', clustername))
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Viruses', clustername))

allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Caudoviricetes', fullquery_clusterLCA))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('phage|Phage', fullquery_clusterLCA))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Carnation_latent', fullquery_clusterLCA))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Microviridae', fullquery_clusterLCA))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Porcine_bastrovirus', fullquery_clusterLCA))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Porcine_picobirnavirus', fullquery_clusterLCA))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('1_XXX', fullquery_clusterLCA))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Viruses', fullquery_clusterLCA))

allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Caudoviricetes', cluster))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('phage|Phage', cluster))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Carnation_latent', cluster))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Microviridae', cluster))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Porcine_bastrovirus', cluster))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Porcine_picobirnavirus', cluster))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('1_XXX', cluster))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Viruses', cluster))

allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Inoviridae', fullquery_clusterLCA))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Leviviridae', fullquery_clusterLCA))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Inoviridae', cluster))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Leviviridae', cluster))
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Cystoviri', fullquery_clusterLCA))
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!grepl('Cystoviri', cluster))


### also check again to see if there are any singleton targets after these filters...
# ## need to run a new group_by() and count
# #allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(clustername) %>% mutate(newcount=n())
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter((grepl('_NODE', query) & clustersize == 1) | clustersize > 1)
# #allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence$newcount <- NULL
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% ungroup()
# 
## different filter - this also will remove clusters of two that are only nttargets - already done before but need to redo for non-phage, now done below
#allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!(bioproject == "NTtarget" & bioprojectsbycluster == 1))


## also get new counts and bioprojects for clusterLCA!!  #allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(clustername) %>% mutate(newcount=n())
## clusterLCA_contigcount & clusterLCA_bioprojectcount
## note save these counts for 'interesting' set

## add in here first a removal of the "NTtarget" bioprojects or convert to NA in a separate dataframe, then use that to summarize a new column clusterLCA_bioprojectcount...
## WE ALSO NEED TO DO THIS FOR bioprojectsbycluster - now done above for full virus set but see below, good to calculate this and clustersize after removing phage
# 
## instead of filtering, convert to NA and use n_distinct(..., na.rm = TRUE) - was a new dataframe allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_notargets but keeping same now and adding replace_na below
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% mutate(bioproject = na_if(bioproject, "NTtarget"))
# ##   bioprojectsbycluster %>% group_by(cluster) %>% mutate(bioprojectsbycluster=n_distinct(bioproject))  
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(cluster) %>% mutate(bioprojectsbycluster2 = n_distinct(bioproject, na.rm = TRUE)) %>% ungroup()
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% relocate(bioprojectsbycluster2, .after = bioprojectsbycluster)

allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(clusterLCA) %>% mutate(clusterLCA_bioprojectcount = n_distinct(bioproject, na.rm = TRUE)) %>% ungroup()
#allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(clusterLCA) %>% mutate(clusterLCA_bioprojectcount=n_distinct(bioproject))

allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(clusterLCA) %>% add_tally() %>% rename(clusterLCA_contigcount = n) %>% ungroup()
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% relocate(clusterLCA_bioprojectcount, .after = clusterLCA)
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% relocate(clusterLCA_contigcount, .after = clusterLCA)

## also calculate a cluster per clusterLCA count
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(clusterLCA) %>% mutate(clusterLCA_clustercount = n_distinct(cluster, na.rm = TRUE)) %>% ungroup()
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% relocate(clusterLCA_clustercount, .after = clusterLCA_bioprojectcount)

## replace old clustername with clusterLCA & count with clustersize_nonphage
## best to use old cluster here not clusterLCA
#allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(clusterLCA) %>% summarize(clustersize_nonphage=n()) %>% arrange(desc(clustersize_nonphage))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(cluster) %>% summarize(clustersize_nonphage=n()) %>% arrange(desc(clustersize_nonphage)) %>% ungroup()
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(clusterLCA) %>% summarize(clusterLCA_contigcount=n()) %>% arrange(desc(clusterLCA_contigcount)) %>% ungroup()
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(clusterLCA) %>% summarize(clusterLCA_clustercount = n_distinct(cluster, na.rm = TRUE))  %>% arrange(desc(clusterLCA_clustercount)) %>% ungroup()

## also here is a good place to recalculate clustersize & bioprojectsbycluster but for non-phage subset...
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(cluster) %>% mutate(bioprojectsbycluster2 = n_distinct(bioproject, na.rm = TRUE)) %>% ungroup()
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(cluster) %>% add_tally() %>% rename(clustersize2 = n) %>% ungroup()

## and filter if clusterLCA_bioprojectcount == 0
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!(clusterLCA_bioprojectcount == 0))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!(clusterLCA_bioprojectcount == 0))

## now add back NTtarget using replace_na() replace_na(list(cpart2 = ""))
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% replace_na(list(bioproject = "NTtarget"))
# ## also remove bioprojectsbycluster & rename bioprojectsbycluster2, same with clustersize & clustersize2
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% select(-bioprojectsbycluster) %>% rename(bioprojectsbycluster = bioprojectsbycluster2)
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% select(-clustersize) %>% rename(clustersize = clustersize2)
## another filter if bioprojectsbycluster == 0
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!(bioprojectsbycluster == 0))

allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% relocate(clustersize, .before = bioproject)
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% relocate(bioprojectsbycluster, .before = bioproject)

## summaries of clusters & clusterLCAs by analysis

# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_LCAsummary <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>%
#   group_by(clusterLCA, analysis_used) %>%
# #  summarise(clustersize_nonphage = n(), .groups = "drop") %>%
#   summarise(clustersize_byLCAandanalysis = n(), clustersize = first(clustersize), clusterLCA_bioprojectcount = first(clusterLCA_bioprojectcount), .groups = "drop") %>%
#   arrange(desc(clustersize),clusterLCA,desc(clustersize_byLCAandanalysis)) %>% select(-clustersize) %>% select(-clusterLCA_bioprojectcount)

## use this but without analysis groupings for treemaps below...
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_LCAsummary <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>%
  group_by(clusterLCA, analysis_used) %>%
  #  summarise(clustersize_nonphage = n(), .groups = "drop") %>%
  summarise(clustersize_byLCAandanalysis = n(), clusterLCA_contigcount = first(clusterLCA_contigcount), clusterLCA_bioprojectcount = first(clusterLCA_bioprojectcount), .groups = "drop") %>%
  arrange(desc(clusterLCA_contigcount),clusterLCA,desc(clustersize_byLCAandanalysis)) %>% ungroup()

allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_summary <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>%
  group_by(cluster, analysis_used) %>%
  #  summarise(clustersize_nonphage = n(), .groups = "drop") %>%
  summarise(clustersize_byanalysis = n(), clusterLCA = first(clusterLCA), clustersize = first(clustersize), bioprojectsbycluster = first(bioprojectsbycluster), .groups = "drop") %>%
  arrange(desc(clustersize),desc(bioprojectsbycluster),clusterLCA, desc(clustersize_byanalysis)) %>% ungroup()

## some naming stats: average length of allchunks_diamondnr_andblastntclustered_viruses_nonphage_withsequence_forfasta$read was
## allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta
## nl <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta %>% mutate(namelength = nchar(read))
# mean(nl$namelength) = 219 median(nl$namelength) = 216
## vs allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence$fullquery_clusterLCA[1]
# nl2 <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% mutate(namelength = nchar(fullquery_clusterLCA)) %>% select(namelength)
# mean(nl2$namelength) = 89.1 median(nl2$namelength) = 81

# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_summary <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>%
#   group_by(clustername, analysis_used) %>%
#   #  summarise(count = n(), .groups = "drop") %>%
#   summarise(count = n(), clustersize = first(clustersize), bioprojectsbycluster = first(bioprojectsbycluster), .groups = "drop") %>%
#   arrange(desc(clustersize),desc(bioprojectsbycluster),clustername, desc(count))



## now at end rename many columns in both allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b and allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence
## cluster to cluster_longestcontigname
## query to contig_name
## shortquery to shortcontigname
## fullquery to contig_withLCA
## fullquery_clusterLCA to contig_withLCA_withcluster

allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% rename(contig_name = query)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% rename(shortcontigname = shortquery)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% rename(contig_withLCA = fullquery)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% rename(contig_withLCA_withcluster = fullquery_clusterLCA)

allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% rename(contig_name = query)
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% rename(shortcontigname = shortquery)
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% rename(contig_withLCA = fullquery)
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% rename(contig_withLCA_withcluster = fullquery_clusterLCA)

## save larger set
## save
write.table(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b, file = paste0("taxonomy_hits_viruses_withsequenceandclusters_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## also always have a generic version saved NOTE THIS WILL CHANGE, BUT IS USED FOR PART 3 AND ALSO VIRUS CURATION...
write.table(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b, file = paste0("taxonomy_hits_viruses_withsequenceandclusters_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## save interesting set
write.table(allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_summary, file = paste0("taxonomy_hits_viruses_nonphage_withsequenceandclusters_analysissummary",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_LCAsummary, file = paste0("taxonomy_hits_viruses_nonphage_withsequenceandLCAclusters_analysissummary",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence, file = paste0("taxonomy_hits_viruses_nonphage_withsequenceandclusters_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
## and most recent
write.table(allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence, file = paste0("taxonomy_hits_viruses_nonphage_withsequenceandclusters_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


###########################
# ## also a treemap but showing cluster info? allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence
#allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_summary

# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(cluster) %>% summarize(clustersize_nonphage=n()) %>% arrange(desc(clustersize_nonphage))
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% group_by(clusterLCA) %>% summarize(clusterLCA_contigcount=n()) %>% arrange(desc(clusterLCA_contigcount))

# viruscluster_heatmap <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>%
#   group_by(clustername) %>%
#   summarise(count = n(), bioprojectsbycluster = first(bioprojectsbycluster)) %>%
#   ungroup()

# removing this treemap only of long cluster names - better to simply show the clusterLCA treemap below - update: adding back but using shortclustername
viruscluster_heatmap0 <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!(bioproject == "NTtarget")) %>% group_by(shortclustername) %>% add_tally() %>% select(-clustersize) %>% rename(clustersize = n) %>% ungroup() %>%
  relocate(clusterLCA_contigcount, .before = clusterLCA_bioprojectcount)


viruscluster_heatmap0 <- viruscluster_heatmap0 %>%
  group_by(shortclustername) %>%
  summarise(clustersize = first(clustersize), bioprojectsbycluster = first(bioprojectsbycluster), .groups = "drop") %>%
  arrange(desc(clustersize),desc(bioprojectsbycluster))



## then sort descending, take just first 25, and use these for treemap
viruscluster_heatmap1 <- viruscluster_heatmap0 %>% arrange(desc(clustersize)) %>% slice_head(n = 25)


viruscluster_heatmap1 <- viruscluster_heatmap1 %>% mutate(shortclustername = str_replace_all(shortclustername, c("Severe acute respiratory syndrome coronavirus 2" = "SARS-CoV-2")))
#viruscluster_heatmap1 <- viruscluster_heatmap1 %>% mutate(shortclustername = str_replace_all(shortclustername, c("Severe_acute_respiratory_syndrome-related_coronavirus" = "SARS-CoV-2")))
viruscluster_heatmap1 <- viruscluster_heatmap1 %>% mutate(shortclustername = str_replace_all(shortclustername, c("Human_immunodeficiency_virus_1" = "HIV-1")))


## for first heatmap, make cluster name easier just for plot - do this by separate at pipe plus a few changes
#viruscluster_heatmap1 <- viruscluster_heatmap1 %>% separate_wider_delim(cluster, delim = "|", names = c("cpart1", "cpart2"), too_many = "drop", too_few = "align_start", cols_remove = FALSE)
## for treemap drop sizes in clustername _delim
viruscluster_heatmap1 <- viruscluster_heatmap1 %>% separate_wider_delim(shortclustername, delim = "|", names = c("shortcluster"), too_many = "drop", too_few = "align_start", cols_remove = FALSE)

# ## new column combining text & numbers, also rewording unite & str_
# ## trim cluster name to 50 chars
# viruscluster_heatmap1 <- viruscluster_heatmap1 %>% mutate(cpart1 = substr(cpart1, 1, 40))
# #viruscluster_heatmap1 <- viruscluster_heatmap1 %>% mutate(cluster = substr(cluster, 1, 50))
# #viruscluster_heatmap1 <- viruscluster_heatmap1 %>% unite(label, clusterLCA, clustersize_nonphage, bioprojectsbycluster, sep = "\n", remove = FALSE)

#library(treemapify)
## updating to get commas in numbers & italics
# viruscluster_heatmap2 <- viruscluster_heatmap2 %>%
#   mutate(
#     formatted_label = paste0(cluster, "\n", (comma(clustersize_nonphage)), "\t", (comma(bioprojectsbycluster)))
#   )
# viruscluster_heatmap1 <- viruscluster_heatmap1 %>% replace_na(list(cpart2 = ""))
# 
# viruscluster_heatmap1 <- viruscluster_heatmap1 %>%
#   mutate(
#     formatted_label = paste0(cpart1, "\n", cpart2, "\n", (comma(clustersize)), "\t", (comma(bioprojectsbycluster)))
#   )

## changing tab below to space - tab breaks in png version...
viruscluster_heatmap1 <- viruscluster_heatmap1 %>%
  mutate(
    formatted_label = paste0(shortcluster, "\n", (comma(clustersize)), " ", (comma(bioprojectsbycluster)))
  )


viruscluster_treemap <- ggplot(viruscluster_heatmap1,
                               aes(area = clustersize, fill = shortcluster, label = formatted_label, subgroup = shortcluster)) +
  geom_treemap(layout = "squarified") +
  geom_treemap_text(place = "centre", size = 18, fontface = "italic") +  # Labels inside blocks
  theme_minimal(base_size = 14, base_family = "sans") +
  labs(title = "Treemap of non-phage viral contig clusters, with contig & bioproject counts per cluster", fill = "Clusters") +
  theme(legend.position = "none")  # Removes the external legend


## UPDATE - REMOVE NTTARGETS FROM THESE (on dataset just before plots)
## ALSO HAVE THIRD VALUE IN THE TREEMAPS - THAT INDICATES THE NUMBER OF CLUSTERS! clusterLCA_clustercount
##############################
## now repeat for clusterLCA!
## first re-calculate clusterLCA_contigcount (other two values are ok)

## remove targetNT, then recalculate clusterLCA_contigcount
#allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% dplyr::filter(!(bioproject == "NTtarget"))

viruscluster_heatmap <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!(bioproject == "NTtarget")) %>% group_by(clusterLCA) %>% add_tally() %>% select(-clusterLCA_contigcount) %>% rename(clusterLCA_contigcount = n) %>% ungroup() %>%
  relocate(clusterLCA_contigcount, .before = clusterLCA_bioprojectcount)

viruscluster_heatmap <- viruscluster_heatmap %>%
  group_by(clusterLCA) %>%
  summarise(clusterLCA_contigcount = first(clusterLCA_contigcount), clusterLCA_bioprojectcount = first(clusterLCA_bioprojectcount), clusterLCA_clustercount = first(clusterLCA_clustercount), .groups = "drop") %>%
  arrange(desc(clusterLCA_contigcount),clusterLCA)

## then sort descending, take just first 25, and use these for treemap
viruscluster_heatmap2 <- viruscluster_heatmap %>% arrange(desc(clusterLCA_contigcount)) %>% slice_head(n = 25)

## new column combining text & numbers, also rewording unite & str_
## trim cluster name to 50 chars
#viruscluster_heatmap2 <- viruscluster_heatmap2 %>% mutate(cluster = substr(cluster, 1, 50))
viruscluster_heatmap2 <- viruscluster_heatmap2 %>% mutate(clusterLCA = str_replace_all(clusterLCA, c("Severe acute respiratory syndrome coronavirus 2" = "SARS-CoV-2")))
#viruscluster_heatmap2 <- viruscluster_heatmap2 %>% mutate(clusterLCA = str_replace_all(clusterLCA, c("Severe acute respiratory syndrome-related coronavirus" = "SARS-CoV-2")))
viruscluster_heatmap2 <- viruscluster_heatmap2 %>% mutate(clusterLCA = str_replace_all(clusterLCA, c("Human immunodeficiency virus 1" = "HIV-1")))

viruscluster_heatmap2 <- viruscluster_heatmap2 %>% unite(label, clusterLCA, clusterLCA_contigcount, clusterLCA_bioprojectcount, clusterLCA_clustercount, sep = "\n", remove = FALSE)

## changing tab below to space - tab breaks in png version...
## updating to get commas in numbers & italics
viruscluster_heatmap2 <- viruscluster_heatmap2 %>%
  mutate(
    formatted_label = paste0(clusterLCA, "\n", (comma(clusterLCA_contigcount)), " ", (comma(clusterLCA_bioprojectcount)), " ", (comma(clusterLCA_clustercount)))
  )

viruscluster_treemapB <- ggplot(viruscluster_heatmap2, 
                               aes(area = clusterLCA_contigcount, fill = clusterLCA, label = formatted_label, subgroup = clusterLCA)) +
  geom_treemap(layout = "squarified") +  
  geom_treemap_text(place = "centre", size = 18, fontface = "italic") +  # Labels inside blocks
  theme_minimal(base_size = 14, base_family = "sans") +
  labs(title = "Treemap of non-phage viral contig LCA clusters, with contig & bioproject & cluster counts per LCA", fill = "LCA Clusters") +
  theme(legend.position = "none")  # Removes the external legend


ggsave(filename = paste("taxonomy_hits_viruses_nonphage_clusters_treemap_",Sys.Date(),".png", sep=""), viruscluster_treemap, width = 18, height = 9, units = "in", limitsize = FALSE)
ggsave(filename = paste("taxonomy_hits_viruses_nonphage_clusters_treemap_",Sys.Date(),".pdf", sep=""), viruscluster_treemap, width = 18, height = 9, units = "in", limitsize = FALSE)

ggsave(filename = paste("taxonomy_hits_viruses_nonphage_LCAclusters_treemap_",Sys.Date(),".png", sep=""), viruscluster_treemapB, width = 18, height = 9, units = "in", limitsize = FALSE)
ggsave(filename = paste("taxonomy_hits_viruses_nonphage_LCAclusters_treemap_",Sys.Date(),".pdf", sep=""), viruscluster_treemapB, width = 18, height = 9, units = "in", limitsize = FALSE)
#########

### making fastas - with new names much simpler, just use contig_withLCA_withcluster
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence

### large set
## for fastas combine query with cluster info
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta$cluster0 <- "cluster"
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta$cluster1 <- "cluster_size"
# #allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta$cluster2 <- "bioprojects_per_cluster"
# ## full query + cluster
# 
# ## full cluster
# #allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta %>% unite(fullcluster, cluster0, clustername, cluster1, clustersize, cluster2, , sep = "|", remove = FALSE)
# 
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta %>% unite(fullqueryandcluster, fullquery, cluster0, clustername, cluster1, clustersize, bioprojectsbycluster, sep = "|", remove = TRUE)
# #allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta$fullqueryandcluster <- gsub('bioprojects_per_cluster\\|','bioprojects_per_cluster_',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta$fullqueryandcluster)
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta$fullqueryandcluster <- gsub('cluster_size\\|','cluster_size_and_bioprojects_',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta$fullqueryandcluster)
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta$fullqueryandcluster <- gsub('cluster\\|','CLUSTER_',allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta$fullqueryandcluster)
# 
# ## then export to fasta
# ## rename fullqueryandcluster to read
# allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta %>% rename(read = fullqueryandcluster)
allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta %>% rename(read = contig_withLCA_withcluster)

writetoFastafaster(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta, paste0("taxonomy_hits_viruses_withsequenceandclusters_",Sys.Date(),".fasta"))
writetoFastafaster(allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence_forfasta, "taxonomy_hits_viruses_withsequenceandclusters_mostrecent.fasta")

### interesting set - with new names much simpler, just use contig_withLCA_withcluster
## EXPORT INTERESTING TO FASTA TOO. allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta$cluster0 <- "cluster"
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta$cluster1 <- "cluster_size"
# #allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta$cluster2 <- "bioprojects_per_cluster"
# ## full query + cluster
# 
# ## full cluster
# #allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta %>% unite(fullcluster, cluster0, clustername, cluster1, clustersize, cluster2, , sep = "|", remove = FALSE)
# 
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta %>% unite(fullqueryandcluster, fullquery, cluster0, clustername, cluster1, clustersize, bioprojectsbycluster, sep = "|", remove = TRUE)
# #allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta$fullqueryandcluster <- gsub('bioprojects_per_cluster\\|','bioprojects_per_cluster_',allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta$fullqueryandcluster)
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta$fullqueryandcluster <- gsub('cluster_size\\|','cluster_size_and_bioprojects_',allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta$fullqueryandcluster)
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta$fullqueryandcluster <- gsub('cluster\\|','CLUSTER_',allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta$fullqueryandcluster)
# 
# ## then export to fasta
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta %>% rename(read = fullqueryandcluster)
allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta %>% rename(read = contig_withLCA_withcluster)
writetoFastafaster(allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta, paste0("taxonomy_hits_viruses_nonphage_withsequenceandclusters_",Sys.Date(),".fasta"))
writetoFastafaster(allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence_forfasta, "taxonomy_hits_viruses_nonphage_withsequenceandclusters_mostrecent.fasta")



## MAKE VERSION OF INTERESTING BUT WITHOUT TARGETS FOR SALMON just filter analysis_used 
### interesting set no targets -  with new names much simpler, just use contig_withLCA_withcluster
allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_andtargets_withsequence %>% dplyr::filter(!analysis_used == "NTtarget")

# allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta$cluster0 <- "cluster"
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta$cluster1 <- "cluster_size"
# 
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta %>% unite(fullqueryandcluster, fullquery, cluster0, clustername, cluster1, clustersize, bioprojectsbycluster, sep = "|", remove = TRUE)
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta$fullqueryandcluster <- gsub('cluster_size\\|','cluster_size_and_bioprojects_',allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta$fullqueryandcluster)
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta$fullqueryandcluster <- gsub('cluster\\|','CLUSTER_',allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta$fullqueryandcluster)
# 
# ## then export to fasta
# allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta %>% rename(read = fullqueryandcluster)
allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence %>% rename(read = contig_withLCA_withcluster)

writetoFastafaster(allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta, paste0("taxonomy_hits_viruses_nonphage_withsequenceandclusters_notargetsforsalmon0_",Sys.Date(),".fasta"))
writetoFastafaster(allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence_forfasta, "taxonomy_hits_viruses_nonphage_withsequenceandclusters_notargetsforsalmon0_mostrecent.fasta")


### make version of all viruses without targets also for salmon
# this file without NTtarget allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b

allchunks_diamondnr_andblastntclustered_viruses_notargets_withsequence2b <- allchunks_diamondnr_andblastntclustered_viruses_andtargets_withsequence2b %>% dplyr::filter(!analysis_used == "NTtarget")

allchunks_diamondnr_andblastntclustered_viruses_notargets_withsequence2b_forfasta <- allchunks_diamondnr_andblastntclustered_viruses_notargets_withsequence2b %>% rename(read = contig_withLCA_withcluster)

writetoFastafaster(allchunks_diamondnr_andblastntclustered_viruses_notargets_withsequence2b_forfasta, paste0("taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_",Sys.Date(),".fasta"))
writetoFastafaster(allchunks_diamondnr_andblastntclustered_viruses_notargets_withsequence2b_forfasta, "taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_mostrecent.fasta")


## also save tsvs for salmon as well - but remove sequence
allchunks_diamondnr_andblastntclustered_viruses_notargets_nosequence2b <- allchunks_diamondnr_andblastntclustered_viruses_notargets_withsequence2b %>% select(-seq.text)
allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_nosequence <- allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_withsequence %>% select(-seq.text)

write.table(allchunks_diamondnr_andblastntclustered_viruses_notargets_nosequence2b, file = paste0("taxonomy_hits_viruses_clusters_notargetsnosequenceforsalmon_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## also always have a generic version saved NOTE THIS WILL CHANGE, BUT IS USED FOR PART 3 AND ALSO VIRUS CURATION...
write.table(allchunks_diamondnr_andblastntclustered_viruses_notargets_nosequence2b, file = paste0("taxonomy_hits_viruses_clusters_notargetsnosequenceforsalmon_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


write.table(allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_nosequence, file = paste0("taxonomy_hits_viruses_nonphage_clusters_notargetsnosequenceforsalmon0_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(allchunks_diamondnr_andblastntclustered_viruses_nonphage_notargets_nosequence, file = paste0("taxonomy_hits_viruses_nonphage_clusters_notargetsnosequenceforsalmon0_mostrecent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


# ## finally save a version with trimmed sequence?
# ## do this in a last script that also includes both bbduk filter and seqkit trimming and ordering for Salmon
# ## TK in last script  - steps to create taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-11_masked_trimmed.fa & taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-11_masked_trimmed1.fa
# fasta_viruses_forsalmon <- read.fasta("taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_2025-09-11_masked_trimmed1.fa")
# ## rename seq.name to fullquery
# fasta_viruses_forsalmon <- fasta_viruses_forsalmon %>% rename(contig_withLCA_withcluster = seq.name)
# allchunks_diamondnr_andblastntclustered_viruses_notargets_withtrimmedsequence <- left_join(allchunks_diamondnr_andblastntclustered_viruses_notargets_nosequence2b, fasta_viruses_forsalmon)
# write.table(allchunks_diamondnr_andblastntclustered_viruses_notargets_withtrimmedsequence, file = paste0("taxonomy_hits_viruses_withsequenceandclusters_notargetsforsalmon_withtrimmedsequence",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


##################################################################

### SUBSETTING FOR CREATING MANY FASTAS, ONE PER CLUSTER is now in separate script...
