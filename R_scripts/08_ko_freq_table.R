# Format DRAM output to create KO frequency tables

# running on bunya server

setwd("/scratch/project/micro_inducers")

library(tidyverse)

## import data ----

# gene annotations
annotations.df <- read.table(
  "/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/08_annotate/dram/dram_annotate_genes/cd_100/annotations.tsv",
  sep = '\t', header = T, strip.white = T, fill = T, quote = ""
)

# gene abundance
genes.df <- read.table(
  "/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/08_annotate/cdhit/rpkm/combined_rpkm_results/all_samples_cdhit100_rpkm_per_gene_normalised.tsv", 
  header = T, sep = "\t", strip.white = T, fill = T
)

# metadata
meta.df <- read.csv(
  "/scratch/project/micro_inducers/data/metadata_files/sample_metadata_filtered.csv",
  header = T, strip.white = T
)

# subset samples to meta
genes.df <- genes.df %>%
  select(Gene_ID, all_of(meta.df$Sample_ID))

# remove any genes no longer in df
genes.df <- genes.df[rowSums(genes.df[,-1]) > 0,]

## create KO table ----

# clean contig ids
names(annotations.df)[1] <- "Gene_ID"
annotations.df$Gene_ID <- gsub("total.protein.faa_", "", annotations.df$Gene_ID)

# select columns
annotations.df <- annotations.df %>%
  select(Gene_ID, ko_id)

# join annotations with coverage
ko.df <- inner_join(
  annotations.df, genes.df, by = "Gene_ID"
)

# remove rows with no ko_id
ko.df <- ko.df %>%
  filter(ko_id != "")

# group by ko
ko.df <- ko.df %>%
  group_by(ko_id) %>%
  summarise(across(-Gene_ID, sum)) %>%
  as.data.frame()

## save
write.table(
  x = ko.df, 
  file ="/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/08_annotate/r/ko_frequency_table.tsv",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t")



