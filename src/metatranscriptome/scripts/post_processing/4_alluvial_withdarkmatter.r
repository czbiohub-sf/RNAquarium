library(data.table)
library(scales)
library(fs)
library(taxonomizr)
library(tidyverse)
library(ggalluvial)
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


## need to load this file before setting directory to outpath
counts <- scan("pipeline_counts.txt")
names(counts) <- c(
  "pre_filter", "post_filter", "nt_hits", 
  "nr_hits", "nt_or_nr", "dark_matter"
)

# Read the integer from count.txt do this at start in main folder!
count_value <- as.numeric(readLines("taxonomy_hits_notfoundinNTorNR_count0.txt"))

outpath <- str_c(workingpath, "/RNAquarium_outputs")
setwd(outpath)

## need to create output folder if it doesn't already exist
dir.create(file.path(outpath,"stats"))
statspath <- str_c(outpath, "/stats")

## change setwd for the two inputs
## first let's load the step 1 alluvial plots

### read in all files first, then set output path to virus subfolder...
alluvialfiles <- fs::dir_ls(glob="taxonomy_hits_nonhost_alluvialplot_counts*", recurse = FALSE)
## chunk names
#alluvialfile <- alluvialfiles[1]
## instead use more recent file
alluvialfile <- alluvialfiles[which.max(file_info(alluvialfiles)$modification_time)]


data_alluvialb0 <- read_tsv(alluvialfile)
#write.table(data_alluviala, file = paste0("taxonomy_hits_nonhost_alluvialplot_counts_",Sys.Date(),".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


#####

#### retain old order with missing last
data_alluviala <- data_alluvialb0 %>%
  mutate(across(everything(), ~replace_na(.x, "Missing")))
#data_alluviala <- data_alluviala %>% rename(NT = taxoncategorysimple_NTclustered)
#data_alluviala <- data_alluviala %>% rename(NR = taxoncategorysimple_NR)
#data_alluviala <- data_alluviala %>% rename(used_NTorNR = taxoncategorysimple_NTorNR)

sumcountonepercent <- (sum(data_alluviala$count) / 100)
## set of order of each category using a factor??
data_alluviala <- data_alluviala %>% arrange(desc(count))

### alternate plan for using treemap order, out to 11
#data_heatmap0 <- data_alluviala %>%
#  group_by(used_NTorNR) %>%
#  summarise(totalcount = sum(count)) %>%
#  ungroup()
#data_heatmap0 <- data_heatmap0 %>% rename(taxoncategorysimple_NTorNR = used_NTorNR)
#data_heatmap0 <- data_heatmap0 %>% rename(count = totalcount)
#alluvialorder11o <- data_heatmap0 %>% arrange(desc(count)) %>% slice_head(n = 11)
#alluvialorder11 <- unique(alluvialorder11o$taxoncategorysimple_NTorNR)
#alluvialorder11[12] <- "Missing"

## better code because sometimes totalcount order isn't order of categories in used_NTorNR
used_category_totals <- data_alluviala %>%
  group_by(used_NTorNR) %>%
  summarise(total_count = sum(count), .groups = 'drop') %>%
  arrange(desc(total_count)) %>%
  slice_head(n = 11)
#top_used_categories <- used_category_totals$used_NTorNR
alluvialorder11 <- used_category_totals$used_NTorNR
alluvialorder11[12] <- "Missing"

## updating to make a full alluvialorder14
used_category_totals <- data_alluviala %>%
  group_by(used_NTorNR) %>%
  summarise(total_count = sum(count), .groups = 'drop') %>%
  arrange(desc(total_count))
#top_used_categories <- used_category_totals$used_NTorNR
alluvialorder14 <- used_category_totals$used_NTorNR
alluvialorder14[15] <- "Missing"

##############################
## then the new dark matter count to add at end - only saving diagram there but also in main script move files to statspath as well
setwd(statspath)

## ADDING IN A SINGLE VALUE FROM A COUNT.TXT FILE...

# Read the integer from count.txt do this at start in main folder!
#count_value <- as.numeric(readLines("taxonomy_hits_notfoundinNTorNR_count0.txt"))

# Create a data frame with one row and the specified column names

# data_alluvialb1 <- data.frame(
#   taxoncategorysimple_NTclustered = "Missing",
#   taxoncategorysimple_NTorNR = "Missing",
#   taxoncategorysimple_NR = "Missing",
#   count = count_value,
#   stringsAsFactors = FALSE
# )



data_alluvialb1 <- data.frame(
  NT = "Missing",
  used_NTorNR = "Missing",
  NR = "Missing",
  count = count_value,
  stringsAsFactors = FALSE
)

# View the resulting data frame
print(data_alluvialb1)


## now alluvial plots & treemaps

setwd(outpath)

data_alluvialadm <- bind_rows(data_alluvialb0,data_alluvialb1)
#rm(data_alluvialb0)
#rm(data_alluvialb1)



#data_alluvialadm <- data_alluvialadm %>% rename(NT = taxoncategorysimple_NTclustered)
#data_alluvialadm <- data_alluvialadm %>% rename(NR = taxoncategorysimple_NR)
#data_alluvialadm <- data_alluvialadm %>% rename(used_NTorNR = taxoncategorysimple_NTorNR)

#sumcountonepercent <- (sum(data_alluviala$count) / 100)
#sumcountonepercentdm <- (sum(data_alluvialadm$count) / 100)
#sumcountonepercent <- (sum(data_alluvialadm$count) / 100)

## set of order of each category using a factor??
data_alluvialadm <- data_alluvialadm %>% arrange(desc(count))


## remove rows below cutoff of sumcountonepercent
sumcountonepercent4 <- sumcountonepercent / 4
data_alluvialdm <- data_alluvialadm %>% dplyr::filter(count > sumcountonepercent4)


## now to figure out 11 hues
#library(colorspace)
cl = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#FF62BC", "#A3A500")
#lightcl <- lighten(cl, amount = 0.3)
#darkcl <- darken(cl, amount = 0.3)
#[1] "#9F1718" "#2D577C" "#397737" "#693870" "#B05601" "#ACAC00" "#743C1E" "#D32092" "#696969" "#CC018C" "#6F7105"
# scale_fill_manual(values = c( "#00969A", "#4D6E00", "#C12600", "#9E00DF")) + ## 0.3
darkcl <- c("#9F1718", "#2D577C", "#397737", "#693870", "#B05601", "#ACAC00", "#743C1E", "#D32092", "#696969", "#CC018C", "#6F7105")
lightcl <- c("#FF6B6B", "#69A4DF", "#6ECD6C", "#C17CCD", "#FFA984", "#FFFF91", "#D4805B", "#FFA7D2", "#B7B7B7", "#FF9ACE", "#BFC144")

## using 12 hues with dark gray last for missing dark matter!

cl12 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#36454F", "#36454F", "#36454F")
darkcl12 <- c("#9F1718", "#2D577C", "#397737", "#693870", "#B05601", "#ACAC00", "#743C1E", "#D32092", "#696969", "#36454F", "#36454F", "#36454F")
lightcl12 <- c("#FF6B6B", "#69A4DF", "#6ECD6C", "#C17CCD", "#FFA984", "#FFFF91", "#D4805B", "#FFA7D2", "#B7B7B7", "#36454F", "#36454F", "#36454F")

## expanded palette with 14 hues - but dark gray last
#lightcl14 <- c("#FF6B6B", "#69A4DF", "#6ECD6C", "#C17CCD", "#FFA984", "#FFFF91", 
#             "#D4805B", "#FFA7D2", "#B7B7B7", "#FF9ACE", "#BFC144",
#             "#B3E5FC", "#FFCC80", "#36454F")  # 3 new colors
cl14 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
             "#A65628", "#F781BF", "#999999", "#FF62BC", "#A3A500",
             "#36454F", "#36454F", "#36454F")  # more gray

#data_alluvialdm$NT <- factor(data_alluvialdm$NT, levels = alluvialorder11, ordered = TRUE)
#data_alluvialdm$NR <- factor(data_alluvialdm$NR, levels = alluvialorder11, ordered = TRUE)
#data_alluvialdm$used_NTorNR <- factor(data_alluvialdm$used_NTorNR, levels = alluvialorder11, ordered = TRUE)
## this works, but also have to remove "decreasing = FALSE" across the commands
## update to 14 levels
data_alluvialdm$NT <- factor(data_alluvialdm$NT, levels = alluvialorder14, ordered = TRUE)
data_alluvialdm$NR <- factor(data_alluvialdm$NR, levels = alluvialorder14, ordered = TRUE)
data_alluvialdm$used_NTorNR <- factor(data_alluvialdm$used_NTorNR, levels = alluvialorder14, ordered = TRUE)



alluvialdm_plotall <- ggplot(data_alluvialdm, aes(axis1 = NT, axis2 = used_NTorNR, axis3 = NR, y = count)) +
  geom_alluvium(aes(fill = used_NTorNR), width = 1/2) +
  stat_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 5, min.y = sumcountonepercent) +
  scale_x_discrete(limits = c("NT", "best (NT or NR)", "NR"), expand = c(0.08, 0.02)) +
  scale_y_continuous(labels = scales::comma) +
  # scale_fill_brewer(palette = "Set1") + ## Set3 has more categories but worse color scheme
  scale_fill_manual(values = cl14) + ## 0.3, instead of darkcl use cl - update changing cl12 to cl14
  theme_minimal(base_size = 14) + ## also theme_minimal, theme_void
  #  theme_minimal(base_family="Helvetica", base_size = 16) + ## also theme_minimal, theme_void
  
  #  labs(title = "Alluvial Diagram of NT vs NR taxonomic categories", fill = "best (NT or NR)",
  labs(title = "Alluvial Diagram of NT vs NR taxonomic categories, including missing Dark Matter", fill = "Taxon",
       x = "", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16))

#alluvialdm_plotall


ggsave(filename = paste("taxonomy_hits_nonhost_alluvialplot_all_withdarkmatter_",Sys.Date(),".png", sep=""), alluvialdm_plotall, width = 18, height = 9, units = "in", limitsize = FALSE)
#ggsave(filename = paste("allchunks_alluvialplot_all0_",Sys.Date(),".png", sep=""), alluvial_plotall, width = 5.6, height = 3.4, units = "in", limitsize = FALSE)
ggsave(filename = paste("taxonomy_hits_nonhost_alluvialplot_all_withdarkmatter_",Sys.Date(),".pdf", sep=""), alluvialdm_plotall, width = 18, height = 9, units = "in", limitsize = FALSE)

################################################
######## final bar plot ##########
used_category_totals <- data_alluviala %>% group_by(used_NTorNR) %>% summarise(n = sum(count)) %>% arrange(desc(n))
## remove mollusca annelida platyhelminthes dplyr::filter(viruscategorysimple_NTorNR != "Phage")
used_category_totals11 <- used_category_totals %>% dplyr::filter(used_NTorNR != "Mollusca")
used_category_totals11 <- used_category_totals11 %>% dplyr::filter(used_NTorNR != "Platyhelminthes")
used_category_totals11 <- used_category_totals11 %>% dplyr::filter(used_NTorNR != "Annelida")

used_category_totals11 <- used_category_totals11 %>%
  mutate(used_NTorNR = factor(used_NTorNR, levels = used_NTorNR))

lightcl11 <- c("#FF6B6B", "#69A4DF", "#6ECD6C", "#C17CCD", "#FFA984", "#FFFF91", 
               "#D4805B", "#FFA7D2", "#FF9ACE", "#BFC144",
               "#B3E5FC")

# ggplot(used_category_totals11, aes(x = used_NTorNR, y = n)) +
#   geom_col() + scale_fill_manual(values = lightcl11) +
#   scale_y_log10(labels = scales::comma) +
#   labs(
#     x = "Category",
#     y = "Count (log scale)",
#     title = "Distribution across taxonomic categories"
#   ) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

contig_barplot <- ggplot(used_category_totals11, aes(x = used_NTorNR, y = n, fill = used_NTorNR)) +
  geom_col() + 
  scale_fill_manual(values = lightcl11) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Taxonomic category",
    y = "Count",
    title = "Distribution of non-host transcripts across taxonomic categories"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")  # Remove legend since x-axis already shows categories

ggsave(filename = paste("taxonomy_hits_nonhost_taxoncategories_barplot_",Sys.Date(),".png", sep=""), contig_barplot, width = 9, height = 9, units = "in", limitsize = FALSE)
ggsave(filename = paste("taxonomy_hits_nonhost_taxoncategories_barplot_",Sys.Date(),".pdf", sep=""), contig_barplot, width = 9, height = 9, units = "in", limitsize = FALSE)

############################################################

### ### ### ### ### ### ### ### ### ### ### ### ### 
### adding little pipeline plot

setwd(statspath)


# Read values from pipeline_counts.txt
#counts <- scan("pipeline_counts.txt")
#names(counts) <- c(
#  "pre_filter", "post_filter", "nt_hits", 
#  "nr_hits", "nt_or_nr", "dark_matter"
#)

# Convert to millions
counts_m <- round(counts / 1e6, 1)

# Function to draw flowchart
draw_flowchart <- function() {
  par(mar = c(1, 1, 1, 1))
  plot(1, type = "n", xlim = c(0, 15), ylim = c(-2, 9), axes = FALSE, xlab = "", ylab = "")

  draw_box <- function(x, y, label, count) {
    rect(x, y, x + 2, y + 1.2, col = "lightgray")
    text(x + 1, y + 0.6, label, cex = 1.1, font = 2)
    text(x + 1, y - 0.2, paste0(count, "M contigs"), cex = 1)
  }

  draw_arrow <- function(x1, y1, x2, y2) {
    arrows(x1, y1, x2, y2, length = 0.12, lwd = 1.5)
  }

  # Boxes and arrows
  draw_box(1, 5, "Assemble\ninto contigs", counts_m["pre_filter"])
  draw_arrow(3, 5.6, 4.5, 5.6)
  draw_box(4.5, 5, "Filter Danio\n& human", counts_m["post_filter"])

  draw_arrow(6.5, 5.6, 8.5, 7)
  draw_box(8.5, 7, "NT blast", counts_m["nt_hits"])

  draw_arrow(6.5, 5.6, 8.5, 3)
  draw_box(8.5, 3, "NR diamond", counts_m["nr_hits"])

  draw_arrow(10.5, 7.6, 12.5, 5.6)
  draw_arrow(10.5, 3.6, 12.5, 5.6)
  draw_box(12.5, 5, "In NT or NR", counts_m["nt_or_nr"])

  draw_box(12.5, 1, "NOT in NT or NR\n(dark matter)", counts_m["dark_matter"])
}

# Output PNG
png("metatranscriptome_pipeline_contignumbers.png", width = 1200, height = 700)
draw_flowchart()
dev.off()

# Output PDF
pdf("metatranscriptome_pipeline_contignumbers.pdf", width = 12, height = 7)
draw_flowchart()
dev.off()
