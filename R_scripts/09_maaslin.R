## MaAsLin2 - Microbiome Multivariable Association with Linear Models

# Use linear models to identify genes associated with settlement

# running on bunya server

library(tidyverse)
library(Maaslin2)

## import data ----

# ko table 
ko.df <- read.table(
  "/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/08_annotate/r/ko_frequency_table.tsv", 
  header = T, sep = "\t", strip.white = T, fill = T
)

# metadata
meta.df <- read.csv(
  "/scratch/project/micro_inducers/data/metadata_files/sample_metadata_filtered.csv",
  header = T, strip.white = T
)

## format data ----

# transpose so sample ID's are rows
rownames(ko.df) <- ko.df$ko_id
ko.df <- ko.df %>%
  select(!ko_id) %>%
  t()

# make sample id row names for meta
rownames(meta.df) <- meta.df$Sample_ID
meta.df <- meta.df %>%
  select(!Sample_ID)

# remove controls 
meta.df <- meta.df %>%
  filter(!Treatment %in% c("control", "negative"))

# create metadata files for each coral
# note: removed samples will automatically be removed from ko df during analysis

# light (1 and 2M) by coral
# psin
meta_psin <- meta.df %>%
  filter(Coral == "P_sinensis")

# dfav
meta_dfav <- meta.df %>%
  filter(Coral == "D_favus")

# easp 
meta_easp <- meta.df %>%
  filter(Coral == "E_aspera") %>%
  filter(Treatment != "D_2M")

# plob
meta_plob <- meta.df %>%
  filter(Coral == "P_lobata") 

# make list
meta.l <- list(
  "psin" = meta_psin,
  "dfav" = meta_dfav,
  "easp" = meta_easp,
  "plob" = meta_plob
)

## run maaslin settlement ----

# sample ids must be row names for meta and feature data

# run list with different models
fit_data <- list()
for(i in seq_along(meta.l)) {
    print(paste0("analysing coral ", names(meta.l)[i]))
    fit_data[[i]] <- Maaslin2(input_data = ko.df,
                              input_metadata = meta.l[[i]],
                              output = paste0("/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/maaslin/", 
                                              names(meta.l)[i], 
                                              "_ko_lm_maaslin"),
                              normalization = "None",
                              transform = "LOG",
                              analysis_method = "LM",
                              random_effects = c("Treatment", "Tank"),
                              fixed_effects = "Percent_settled",
                              min_abundance = 0.001,
                              min_prevalence = 0.1,
                              plot_heatmap = TRUE,
                              cores = 20)
}

## run maaslin treatment ----

# compare treatments to L2M
fit_data <- Maaslin2(input_data = ko.df,
                     input_metadata = meta.df,
                     output = paste0("/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/maaslin/treatment_ko_lm_maaslin"),
                     normalization = "None",
                     transform = "LOG",
                     analysis_method = "LM",
                     random_effects = c("Tank", "Coral"),
                     fixed_effects = "Treatment",
                     reference = c("Treatment,L_2M"),
                     min_abundance = 0.001,
                     min_prevalence = 0.1,
                     plot_heatmap = TRUE,
                     cores = 20)







