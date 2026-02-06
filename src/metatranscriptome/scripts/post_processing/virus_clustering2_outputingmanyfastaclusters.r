library(data.table)
library(scales)
library(fs)
library(taxonomizr)
library(tidyverse)

library(phylotools)
#library(treemapify)


###
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

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

workingpath <- getwd()
workingpathdash <- str_c(workingpath, "/")
diamondpath <- str_c(workingpath, "/nr_diamond")
blastpath <- str_c(workingpath, "/nt_blast")

## saving all outputs to RNAquarium_outputs/virus_outputs
outpath <- str_c(workingpath, "/RNAquarium_outputs")
#setwd(outpath)

outpathvirus <- str_c(outpath, "/virus_outputs")
setwd(outpathvirus)

## for clustering, we want to create new cluster directories for final steps
#clusters_dir_name <- paste0("clusters_", format(Sys.Date(), "%m%d%y"))
clusters_dir_name <- paste0("clusters_forminimap")
dir.create(clusters_dir_name)
outpath2 <- str_c(outpathvirus, "/clusters_forminimap")

## second cluster folder just for phage (not a subfolder)
phage_dir_name <- paste0("clusters_forminimap_phage")
dir.create(phage_dir_name)
outpath2phage <- str_c(outpathvirus, "/clusters_forminimap_phage")

####################################
## after mmseqs2 easy-cluster R commands to curate the results


# clusters_viruses_fasta <- read_tsv("taxonomy_hits_viruses_withsequenceandclusters_mostrecent.tsv")
# 
# clusters_viruses_fasta <- clusters_viruses_fasta %>% select(contig_withLCA_withcluster,seq.text,shortclustername,clustersize)
# clusters_viruses_fasta <- clusters_viruses_fasta %>% rename(read = contig_withLCA_withcluster) %>% rename(cluster = shortclustername) %>% rename(n = clustersize)
# 
# clusters_viruses_fasta_smaller <- clusters_viruses_fasta %>% dplyr::filter(n > 9)
# 
# clusters_viruses_fasta_smaller2 <- clusters_viruses_fasta %>% dplyr::filter(n < 10)
# clusters_viruses_fasta_smaller2 <- clusters_viruses_fasta_smaller2 %>% dplyr::filter(n > 1)
# 
# clusters_viruses_fasta_singletons <- clusters_viruses_fasta %>% dplyr::filter(n == 1)
# 
# 
# ### some summaries then save for combining with allchunks_blastnanddiamond_hits_viruses_withsequence_mostrecent.tsv
# n_distinct(clusters_viruses_fasta_smaller$cluster)
# n_distinct(clusters_viruses_fasta_smaller2$cluster)
# n_distinct(clusters_viruses_fasta_singletons$cluster)

## new split - for all non-singletons, divide into phage & non-phage, save phage in subfolder!

clusters_viruses_fasta <- read_tsv("taxonomy_hits_viruses_withsequenceandclusters_mostrecent.tsv")

clusters_viruses_fasta <- clusters_viruses_fasta %>% select(contig_withLCA_withcluster,seq.text,shortclustername,clustersize,viruscategorysimple_NTorNR,cluster,clusterLCA)
clusters_viruses_fasta <- clusters_viruses_fasta %>% rename(cluster0 = cluster)
#clusters_viruses_fasta <- clusters_viruses_fasta %>% select(contig_withLCA_withcluster,seq.text,shortclustername,clustersize,viruscategorysimple_NTorNR)
clusters_viruses_fasta <- clusters_viruses_fasta %>% rename(read = contig_withLCA_withcluster) %>% rename(cluster = shortclustername) %>% rename(n = clustersize)

clusters_viruses_fasta_singletons <- clusters_viruses_fasta %>% dplyr::filter(n == 1)
clusters_viruses_fasta <- clusters_viruses_fasta %>% dplyr::filter(n > 1)

## smaller to non-phage
## smaller2 to phage
# clusters_viruses_fasta_smaller2 <- clusters_viruses_fasta %>% dplyr::filter(viruscategorysimple_NTorNR != "Non-phage")
# clusters_viruses_fasta_smaller <- clusters_viruses_fasta %>% dplyr::filter(viruscategorysimple_NTorNR == "Non-phage")

# instead of splitting this way, going back to earlier steps to classify these into "Phage" - but we need to include because this script starts with full set, and last greps are not in categorysimple
df_split <- clusters_viruses_fasta %>%
  mutate(group = if_else(
    viruscategorysimple_NTorNR != "Non-phage" | grepl("Caudoviricetes|Prokaryotic", cluster0) | grepl('Caudoviricetes|Prokaryotic', clusterLCA) | grepl('phage|Phage', cluster0) | grepl('phage|Phage', clusterLCA) | grepl('Levivir|levivir|Porcine picobirnavirus', cluster0) | grepl('Levivir|levivir|Porcine picobirnavirus', cluster) | grepl('Inoviridae|Cystoviridae', cluster0) | grepl('Inoviridae|Cystoviridae', clusterLCA),  # Define df_no condition
    "phage",
    "nonphage"
  ))

clusters_viruses_fasta_smaller2 <- df_split %>% filter(group == "phage") %>% select(-group)  # Everything else
clusters_viruses_fasta_smaller  <- df_split %>% filter(group == "nonphage") %>% select(-group)   # Matches the 'no' condition


clusters_viruses_fasta_singletons <- clusters_viruses_fasta_singletons %>% select(read,seq.text,cluster,n)
clusters_viruses_fasta_smaller <- clusters_viruses_fasta_smaller %>% select(read,seq.text,cluster,n)
clusters_viruses_fasta_smaller2 <- clusters_viruses_fasta_smaller2 %>% select(read,seq.text,cluster,n)

n_distinct(clusters_viruses_fasta_smaller$cluster)
n_distinct(clusters_viruses_fasta_smaller2$cluster)
n_distinct(clusters_viruses_fasta_singletons$cluster)


## summary steps are in first clustering script - this script goes straight to creating fastas

##################################################################

### BELOW ARE SUBSETTING FOR CREATING MANY FASTAS, ONE PER CLUSTER
## AT THIS POINT USE NEW OUTPUT FOLDER TO SAVE ALL OF THESE...
#setwd(clusters_dir_name)
setwd(outpath2)

## for first batch use new folder - save all in a single folder, makes next minimap2 step easier
## see better solution above
#outpath1 <- str_c(workingpath, "/RNAquarium_outputs/virusclusters")
#setwd(outpath1)


# clusters_viruses_smaller <- inner_join(clusters_viruses, clusters_viruses_tally)
# ## now rename & another inner join
# clusters_viruses_fasta_smaller <- inner_join(clusters_viruses_fasta,clusters_viruses_smaller)

## also remove wonky characters (, ), \ replace all with dashes - no to underscores!
## from both cluster & read gsub
clusters_viruses_fasta_smaller$cluster <- gsub("\\(", "_", clusters_viruses_fasta_smaller$cluster)
clusters_viruses_fasta_smaller$cluster <- gsub("\\)", "_", clusters_viruses_fasta_smaller$cluster)
clusters_viruses_fasta_smaller$cluster <- gsub("\\\\", "_", clusters_viruses_fasta_smaller$cluster)
clusters_viruses_fasta_smaller$cluster <- gsub("/", "_", clusters_viruses_fasta_smaller$cluster)

clusters_viruses_fasta_smaller$read <- gsub("\\(", "_", clusters_viruses_fasta_smaller$read)
clusters_viruses_fasta_smaller$read <- gsub("\\)", "_", clusters_viruses_fasta_smaller$read)
clusters_viruses_fasta_smaller$read <- gsub("\\\\", "_", clusters_viruses_fasta_smaller$read)
clusters_viruses_fasta_smaller$read <- gsub("/", "_", clusters_viruses_fasta_smaller$read)
## need to shorten only very long names in cluster - greater than 240, also in cluster
clusters_viruses_fasta_smaller <- clusters_viruses_fasta_smaller %>% mutate(read = substr(read, 1, 240))
clusters_viruses_fasta_smaller <- clusters_viruses_fasta_smaller %>% mutate(cluster = substr(cluster, 1, 240))
## now after trimming long names, remove trailing _ or trailing |NA. grepl("pattern$", taxname_lca_NTorNR):
#clusters_viruses_fasta_smaller$read <- gsub('\\|NA$','',clusters_viruses_fasta_smaller$read)
clusters_viruses_fasta_smaller$read <- gsub('NA$','',clusters_viruses_fasta_smaller$read)
clusters_viruses_fasta_smaller$cluster <- gsub('NA$','',clusters_viruses_fasta_smaller$cluster)
clusters_viruses_fasta_smaller$read <- gsub('_$','',clusters_viruses_fasta_smaller$read)
clusters_viruses_fasta_smaller$cluster <- gsub('_$','',clusters_viruses_fasta_smaller$cluster)

## also can arrange by length??
clusters_viruses_fasta_smaller$length <- nchar(clusters_viruses_fasta_smaller$seq.text)
## final arrange
clusters_viruses_fasta_smaller <- clusters_viruses_fasta_smaller %>% arrange(desc(n),cluster,desc(length)) %>% group_by(cluster)
clusters_viruses_fasta_smaller$name_length <- nchar(clusters_viruses_fasta_smaller$read)
clusters_viruses_fasta_smaller$cluster_length <- nchar(clusters_viruses_fasta_smaller$cluster)
# clusters_viruses_fasta_smaller <- clusters_viruses_fasta_smaller %>% dplyr::filter(name_length >= 240)
# clusters_viruses_fasta_smaller <- clusters_viruses_fasta_smaller %>% dplyr::filter(n >= 2500)


## we want to make cluster a factor, descending by n then name
clusterorder <- clusters_viruses_fasta_smaller %>% arrange(desc(n),cluster,desc(length))
clusterorder2 <- unique(clusterorder$cluster)
#clusterorder <- clusterorder[["cluster"]]

clusters_viruses_fasta_smaller$cluster <- factor(clusters_viruses_fasta_smaller$cluster, levels = clusterorder2, ordered = TRUE)
#clusters_viruses_fasta_smaller %>% group_by(cluster) %>% slice_head(n = 1) %>% select(cluster,n,length)

## will have to shorten very long names

clusters_viruses_fasta_smaller %>%
  group_by(cluster) %>% 
  group_split() %>%
  setNames(unique(clusters_viruses_fasta_smaller$cluster)) %>% 
  imap(~writetoFastafaster(.x, paste0(.y, '.fasta')))


######### REPEAT FOR 2-9 clusters

# 
# clusters_viruses_smaller2 <- inner_join(clusters_viruses, clusters_viruses_tally2)
# ## now rename & another inner join
# clusters_viruses_fasta_smaller2 <- inner_join(clusters_viruses_fasta,clusters_viruses_smaller2)

## also remove wonky characters (, ), \ replace all with dashes - no to underscores!
## from both cluster & read gsub
clusters_viruses_fasta_smaller2$cluster <- gsub("\\(", "_", clusters_viruses_fasta_smaller2$cluster)
clusters_viruses_fasta_smaller2$cluster <- gsub("\\)", "_", clusters_viruses_fasta_smaller2$cluster)
clusters_viruses_fasta_smaller2$cluster <- gsub("\\\\", "_", clusters_viruses_fasta_smaller2$cluster)
clusters_viruses_fasta_smaller2$cluster <- gsub("/", "_", clusters_viruses_fasta_smaller2$cluster)

clusters_viruses_fasta_smaller2$read <- gsub("\\(", "_", clusters_viruses_fasta_smaller2$read)
clusters_viruses_fasta_smaller2$read <- gsub("\\)", "_", clusters_viruses_fasta_smaller2$read)
clusters_viruses_fasta_smaller2$read <- gsub("\\\\", "_", clusters_viruses_fasta_smaller2$read)
clusters_viruses_fasta_smaller2$read <- gsub("/", "_", clusters_viruses_fasta_smaller2$read)
## need to shorten only very long names in cluster - greater than 240, also in cluster
clusters_viruses_fasta_smaller2 <- clusters_viruses_fasta_smaller2 %>% mutate(read = substr(read, 1, 240))
clusters_viruses_fasta_smaller2 <- clusters_viruses_fasta_smaller2 %>% mutate(cluster = substr(cluster, 1, 240))
## now after trimming long names, remove trailing _ or trailing |NA. grepl("pattern$", taxname_lca_NTorNR):
#clusters_viruses_fasta_smaller2$read <- gsub('\\|NA$','',clusters_viruses_fasta_smaller2$read)
clusters_viruses_fasta_smaller2$read <- gsub('NA$','',clusters_viruses_fasta_smaller2$read)
clusters_viruses_fasta_smaller2$cluster <- gsub('NA$','',clusters_viruses_fasta_smaller2$cluster)
clusters_viruses_fasta_smaller2$read <- gsub('_$','',clusters_viruses_fasta_smaller2$read)
clusters_viruses_fasta_smaller2$cluster <- gsub('_$','',clusters_viruses_fasta_smaller2$cluster)

## also can arrange by length??
clusters_viruses_fasta_smaller2$length <- nchar(clusters_viruses_fasta_smaller2$seq.text)
## final arrange
clusters_viruses_fasta_smaller2 <- clusters_viruses_fasta_smaller2 %>% arrange(desc(n),cluster,desc(length)) %>% group_by(cluster)
clusters_viruses_fasta_smaller2$name_length <- nchar(clusters_viruses_fasta_smaller2$read)
clusters_viruses_fasta_smaller2$cluster_length <- nchar(clusters_viruses_fasta_smaller2$cluster)

## we want to make cluster a factor, descending by n then name
clusterorder <- clusters_viruses_fasta_smaller2 %>% arrange(desc(n),cluster,desc(length))
clusterorder2 <- unique(clusterorder$cluster)
#clusterorder <- clusterorder[["cluster"]]

clusters_viruses_fasta_smaller2$cluster <- factor(clusters_viruses_fasta_smaller2$cluster, levels = clusterorder2, ordered = TRUE)


### for phage - use separate subfolder!
setwd(outpath2phage)

clusters_viruses_fasta_smaller2 %>%
  group_by(cluster) %>% 
  group_split() %>%
  setNames(unique(clusters_viruses_fasta_smaller2$cluster)) %>% 
  imap(~writetoFastafaster(.x, paste0(.y, '.fasta')))


## finally repeat for singletons
# clusters_viruses_singletons <- inner_join(clusters_viruses, clusters_viruses_singletons)
# ## now rename & another inner join
# clusters_viruses_fasta_singletons <- inner_join(clusters_viruses_fasta,clusters_viruses_singletons)

## also remove wonky characters (, ), \ replace all with dashes - no to underscores!
## from both cluster & read gsub
clusters_viruses_fasta_singletons$cluster <- gsub("\\(", "_", clusters_viruses_fasta_singletons$cluster)
clusters_viruses_fasta_singletons$cluster <- gsub("\\)", "_", clusters_viruses_fasta_singletons$cluster)
clusters_viruses_fasta_singletons$cluster <- gsub("\\\\", "_", clusters_viruses_fasta_singletons$cluster)
clusters_viruses_fasta_singletons$cluster <- gsub("/", "_", clusters_viruses_fasta_singletons$cluster)

clusters_viruses_fasta_singletons$read <- gsub("\\(", "_", clusters_viruses_fasta_singletons$read)
clusters_viruses_fasta_singletons$read <- gsub("\\)", "_", clusters_viruses_fasta_singletons$read)
clusters_viruses_fasta_singletons$read <- gsub("\\\\", "_", clusters_viruses_fasta_singletons$read)
clusters_viruses_fasta_singletons$read <- gsub("/", "_", clusters_viruses_fasta_singletons$read)



## first remove the targets from list of singletons... - filter to only include _NODE  %>% dplyr::filter(grepl('_NODE', read)) OR S_NODE
#clusters_viruses_fasta_singletons <- clusters_viruses_fasta_singletons %>% dplyr::filter(grepl('_NODE_', read))

## just save as a large fasta...
## for the one singleton fasta - change working directory:
setwd(outpathvirus)
#writetoFastafaster(clusters_viruses_fasta_singletons,"RNaquarium_allchunks_blastnanddiamond_hits_viruses_interestinglist_withsequence2024-08-27_dedupedrenamed_singletons.fasta")
writetoFastafaster(clusters_viruses_fasta_singletons, paste0("taxonomy_hits_viruses_clustered_singletons_",Sys.Date(),".fasta"))

## then finally minimap2 will be yet separate...
#module load minimap2/2.26

