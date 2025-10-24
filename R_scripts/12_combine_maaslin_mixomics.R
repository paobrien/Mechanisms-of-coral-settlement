## identify consistent genes across maaslin and mixomics results

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(KEGGREST)

## import data ----

# ko frequency table
ko.df <- read.table(
  "Data/Metagenomics_full/Raw/dram/genes/ko_frequency_table.tsv",
  sep = "\t", header = T, strip.white = T
)

# metadata
meta.df <- read.csv(
  "Data/Metagenomics_full/Formatted/sample_metadata_filtered.csv",
  header = T, strip.white = T
)

# ko annotations
annotations.df <- read.table(
  "/Users/uqpobri2/Desktop/UQ_postdoc/Bioinformatics/MG_biofilm/08_annotate/dram_annotate_genes/ko_annotations.tsv", 
  sep = "\t", header = T, fill = TRUE, quote = "", strip.white = T
)

## import maaslin results
# file path
fp <- "/Users/uqpobri2/Desktop/UQ_postdoc/Bioinformatics/MG_biofilm/10_maaslin/genes/no_ctrl/" 

# directory list
dirs <- list.files(fp, pattern = "_maaslin")

# import significant results from each folder
maas_df.l <- list()
for (i in seq_along(dirs)) {
  maas_df.l[[i]] <- read.table(paste0(fp, dirs[i], "/significant_results.tsv"),
                               header = T,
                               sep = "\t")
  names(maas_df.l)[i] <- gsub("_ko_lm_maaslin", "", dirs[i])
}

## import mixomics results
# file path
fp <- "Results/metagenomics/mixOmics/pls/"

# file list
flist <- list.files(fp, pattern = "spls1_selected_variables")

# import
mix.l <- list()
for (i in seq_along(flist)) {
  mix.l[[i]] <- readRDS(paste0(fp, flist[i]))
  names(mix.l)[i] <- gsub(".rds", "", flist[i]) 
}

## format data ----

# remove contig id
annotations.df <- annotations.df %>%
  select(!contig_id)

# remove duplicate rows
annotations.df <- annotations.df[!(duplicated(annotations.df)),]

# remove blanks and controls
meta.df <- meta.df %>%
  filter(Treatment != "negative") %>%
  filter(Treatment != "control")

ko.df <- ko.df %>%
  select(ko_id, meta.df$Sample_ID)

# combine components from mixomics
comp1_df.l <- list()
comp2_df.l <- list()
for (i in seq_along(mix.l$spls1_selected_variables_c1)) {
  comp1_df.l[[i]] <- data.frame(
    KO = mix.l$spls1_selected_variables_c1[[i]], 
    Comp = rep("Comp1", length(mix.l$spls1_selected_variables_c1[[i]]))
    )
  names(comp1_df.l)[i] <- names(mix.l$spls1_selected_variables_c1)[i]
  comp2_df.l[[i]] <- data.frame(
    KO = mix.l$spls1_selected_variables_c2[[i]], 
    Comp = rep("Comp2", length(mix.l$spls1_selected_variables_c2[[i]]))
    )
  names(comp2_df.l)[i] <- names(mix.l$spls1_selected_variables_c2)[i]
}

# check for duplicated KOs
lapply(comp1_df.l, nrow)
lapply(comp1_df.l, function(x) {
  length(unique(x$KO))
})

## compare results ----

# check which genes are significantly associated with settlement across both analyses
# using component 1 since this had the best separation by settlement
common_genes.l <- list()
for (i in seq_along(comp1_df.l)) {
  common_genes.l[[i]] <- intersect(maas_df.l[[i]]$feature, comp1_df.l[[i]]$KO) 
  names(common_genes.l)[i] <- names(comp1_df.l)[i]
}
lapply(common_genes.l, length)

# repeat but with top n genes from maaslin (thousands of genes, qval not useful cutoff as very different for each coral)
# get top 500 genes
maas_top.l <- lapply(maas_df.l, function(x) {
  x %>% slice_head(n = 500)
})

# combine with above
common_genes.l <- list()
for (i in seq_along(maas_top.l)) {
  common_genes.l[[i]] <- intersect(maas_top.l[[i]]$feature, comp1_df.l[[i]]$KO)
  names(common_genes.l)[i] <- names(comp1_df.l)[i]
}
lapply(common_genes.l, length)

# approx half mixomics genes shared in maaslin analysis. Except Psin where most genes were shared 

# check if any genes consistent among coral
Reduce(intersect, common_genes.l) # none - shared genes across coral occur in maaslin analysis only

# checking mixomics only
columns_list <- lapply(comp1_df.l, function(df) df$KO)
Reduce(intersect, columns_list) 

# checking maaslin only
columns_list <- lapply(maas_df.l, function(df) df$feature)
Reduce(intersect, columns_list) 


## subset to common genes ----

# subset data frames
ko.l <- list()
for (i in seq_along(common_genes.l)) {
  ko.l[[i]] <- ko.df %>%
    filter(ko_id %in% common_genes.l[[i]])
  names(ko.l)[i] <- names(common_genes.l)[i]
}
lapply(ko.l, dim)

# subset samples
# first subset metadata
meta_d <- meta.df %>%
  filter(Coral == "D_favus")

meta_e <- meta.df %>%
  filter(Coral == "E_aspera") %>%
  filter(Treatment != "D_2M")

meta_pl <- meta.df %>%
  filter(Coral == "P_lobata")

meta_ps <- meta.df %>%
  filter(Coral == "P_sinensis")

# combine to list
# metadata
meta.l <- list(
  "meta_dfav" = meta_d,
  "meta_easp" = meta_e,
  "meta_plob" = meta_pl,
  "meta_psin" = meta_ps
)

# select samples
for (i in seq_along(ko.l)) {
  ko.l[[i]] <- ko.l[[i]] %>%
    select(ko_id, meta.l[[i]]$Sample_ID)
}

# make ko id row names
for (i in seq_along(ko.l)) {
  rownames(ko.l[[i]]) <- ko.l[[i]]$ko_id
  ko.l[[i]] <- ko.l[[i]] %>%
    select(!ko_id)
}


## get annotations ----

# combine KOs with multiple annotations
annotations.df <- annotations.df %>%
  group_by(ko_id) %>%
  summarize(kegg_hit = paste(kegg_hit, collapse = "; "))

# remove E ids
annotations.df$kegg_hit <- gsub("\\[.*?\\]", "", annotations.df$kegg_hit)

# remove quotations
annotations.df$kegg_hit <- gsub("[\"']", "", annotations.df$kegg_hit)

# df to plot with annotations
ko_an.l <- ko.l
for (i in seq_along(ko_an.l)) {
  # get ko column
  ko_an.l[[i]]$ko_id <- rownames(ko_an.l[[i]])
  # join with annotation
  ko_an.l[[i]] <- left_join(ko_an.l[[i]], annotations.df, by = "ko_id")
  # unite kegg columns
  ko_an.l[[i]] <- ko_an.l[[i]] %>%
    unite("annotation", ko_id:kegg_hit, sep = ": ", remove = TRUE)
  # change row names
  rownames(ko_an.l[[i]]) <- ko_an.l[[i]]$annotation
  # remove column
  ko_an.l[[i]] <- ko_an.l[[i]] %>%
    select(!annotation)
}



## plot heatmap ----

# for plob heatmap
# remove rows with NA annotation
plob_t <- ko_an.l$ko_plob
plob_t <- plob_t[!grepl("NA", rownames(plob_t)), ]

# log transform
plob_t <- log10(plob_t + 0.1)

# edit long row names
rownames(plob_t)[rownames(plob_t) == "K03532: trimethylamine-N-oxide reductase (cytochrome c), cytochrome c-type subunit TorC"] <- "K03532: trimethylamine-N-oxide reductase (cytochrome c)"
rownames(ko_an.l$ko_easp)[rownames(ko_an.l$ko_easp) == "K11752: diaminohydroxyphosphoribosylaminopyrimidine deaminase / 5-amino-6-(5-phosphoribosylamino)uracil reductase "] <- "K11752: diaminohydroxyphosphoribosylaminopyrimidine deaminase"

# format metadata
for(i in seq_along(meta.l)) {
  rownames(meta.l[[i]]) <- meta.l[[i]]$Sample_ID
}

# get metadata labels
for(i in seq_along(meta.l)) {
  meta.l[[i]]$Treatment <- factor(meta.l[[i]]$Treatment,
                           levels = c("D_2M", "L_1M", "L_2M"),
                           labels = c("Dark 2M", "Light 1M", "Light 2M"))
}

# select variables
meta_ha.l <- lapply(meta.l, function(df) {
  df <- df %>%
    select(Treatment, Percent_settled)
})

# settlement cont
clrs <- viridis(8, option = "plasma")
show_col(clrs)
s_clrs <- colorRamp2(c(0, 50, 100), c("#0D0887FF", "#B93289FF", "#F0F921FF"))

# get column annotations
ha <- HeatmapAnnotation(
  df = meta_ha.l$meta_dfav,
  col = list(
    Treatment = c("Light 2M" = "#FABA39FF", "Light 1M" = "#1AE4B6FF"),
    #Settlement = c("High" = "#000004FF", "Low" = "#B63679FF", "None" = "#FCFDBFFF"),
    Percent_settled = s_clrs)
)

# plot
heat_cols <- viridis(10, option = "turbo")
ht <- Heatmap(ko_an.l$ko_dfav,
              name = "nRPKM (log)",
              col = heat_cols,
              top_annotation = ha, 
              show_column_names = FALSE,
              row_names_max_width = unit(21, "cm"),
              column_split = factor(as.character(meta.l$meta_dfav$Settlement), 
                                    levels = c("High", "None", "Low")),
              column_dend_reorder = TRUE,
              row_km = 2,
              row_km_repeats = 50,
              border = TRUE)
draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

# to save
svg("dfav_settlement_genes_combined.svg", width = 15, height = 10)

draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

dev.off()










