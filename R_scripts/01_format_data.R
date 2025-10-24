## Format metadata, checkm, coverm and gtdb files

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)

## import data ----

# GTDB classification
tax_b <- read.table("Data/Metagenomics_full/Raw/gtdb/derep_bins.bac120.summary.tsv", 
                    sep = '\t', header = T, strip.white = T)

tax_a <- read.table("Data/Metagenomics_full/Raw/gtdb/derep_bins.ar53.summary.tsv", 
                    sep = '\t', header = T, strip.white = T)

tax <- rbind(tax_a, tax_b)

# checkm results
checkm <- read.table("Data/Metagenomics_full/Raw/checkm/derep_bins_checkm.tsv", 
                     sep = '\t', header = T, strip.white = T)

# coverm results
coverm <- read.table("Data/Metagenomics_full/Raw/coverm/dereplicated_genomes_coverage.tsv", 
                     sep = '\t', header = T, strip.white = T)

# metadata
meta.df <- read.csv("Data/Metagenomics_full/Formatted/sample_metadata_all.csv",
                 header = T, strip.white = T)

## format data and make mag metadata file ----

# format taxonomy
tax <- separate(tax, classification, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                sep = ";", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn")

# format genome size
checkm$Genome_size_Mbp <- checkm$Genome_size_bp / 1000000

# combine GTDB with checkm
full.df <- full_join(tax, checkm, join_by("user_genome" == "Bin_Id"))

# add coverm results
full.df <- full_join(full.df, coverm, join_by("user_genome" == "Genome"))
full.df <- full.df %>%
  rename(MAG_ID = user_genome)

# subset to relevant columns
sub.df <- full.df %>% 
  select(MAG_ID, Domain, Phylum, Class, Order, Family, Genus, 
         Species, Completeness, Contamination, Genome_size_Mbp, 
         GC, starts_with("SE"))

## remove poorly sequenced samples (very little data, no mapping)
to_rm <- c("SE2466", "SE2467", "SE2471", "SE2472", "SE2483")

# metadata
meta.df <- meta.df %>%
  filter(!meta.df$Sample_ID %in% to_rm)

# mag data
sub.df <- sub.df %>%
  select(!all_of(to_rm))

# remove data from rubble study
# metadata
rubble <- meta.df %>%
  filter(Tank == "rubble")

meta.df <- meta.df %>%
  filter(Tank != "rubble")

# mag data
sub.df <- sub.df %>%
  select(!all_of(rubble$Sample_ID))

# remove mags no longer in df
# get mags with 0% abundance
sub.df$total <- rowSums(sub.df[ , grepl("^SE", names(sub.df)) ], na.rm = TRUE)

# remove from df
sub.df <- sub.df %>%
  filter(total > 0) %>%
  select(!total)

# clean column names
colnames(sub.df) <- gsub("c000", "", colnames(sub.df))

## save dataframes ----

# summary table
write.table(sub.df, "Data/Metagenomics_full/Formatted/dereplicated_mags_summary_data_filtered.tsv", 
            sep = "\t", row.names = F)

# metadata file
write.csv(meta_f, "Data/Metagenomics_full/Formatted/sample_metadata_filtered.csv", 
          row.names = F)






















