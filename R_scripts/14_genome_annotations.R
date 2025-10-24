# Identify which MAGs have genes of interest

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)

## import data ----

# mags
mags.df <- read.table(
  "Data/Metagenomics_full/Formatted/dereplicated_mags_summary_data_filtered.tsv",
  sep = '\t', header = T, strip.white = T
)

# metadata
meta.df <- read.csv(
  "Data/Metagenomics_full/Formatted/sample_metadata_filtered.csv",
  header = T, strip.white = T
)

# annotations
annotations.df <- read.table(
  "Data/Metagenomics_full/Raw/dram/genomes/annotations_sparse.tsv",
  header = T, strip.white = T, sep = "\t", fill = TRUE
)

# annotations distilled (dram)
distilled.df <- read.table(
  "Data/Metagenomics_full/Raw/dram/genomes/product.tsv", 
  header = T, strip.white = T, sep = "\t", check.names = F
)

# genes of interest
genes.df <- read.csv(
  "Data/Metagenomics_full/Formatted/genes_of_interest.csv",
  header = T, strip.white = T
)

# blast results
tbp <- read.table(
  "Data/Metagenomics_full/Raw/blast/tblastx_significant_results_tbp_to_mags_e06-100.txt", 
  header = T, 
  sep = "\t"
)

cyclo <- read.table(
  "Data/Metagenomics_full/Raw/blast/tblastx_significant_results_cyclo_full_to_mags_e06-100.txt", 
  header = T, 
  sep = "\t"
)

## format data ----

## mags

# remove unmapped data
mags.df <- mags.df %>%
  filter(MAG_ID != "unmapped")

# get taxonomy
tax.df <- mags.df %>% 
  select(MAG_ID, Domain, Phylum, Class, Order, Family, Genus, Species)

# subset to sample rel abundance
mags.df <- mags.df %>%
  select(MAG_ID, starts_with("SE")) 

# clean colnames
colnames(mags.df) <- gsub("c000", "", colnames(mags.df))

## genes

# get number of genes in pathway
genes.df <- genes.df %>%
  group_by(Pathway) %>%
  mutate(Total_genes = n()) %>%
  ungroup()

# update pathways (some genes redundant - hence total genes needed overstated)
genes.df <- genes.df %>%
  mutate(
    Total_genes = case_when(
      Pathway == "Nitrogen fixation" ~ 4,
      Pathway == "Carotenoid (Beta-carotene pathway)" ~ 6,
      Pathway == "Nitric oxide synthesis" ~ 1,
      TRUE ~ Total_genes
    )
  )

## annotations

# subset to MAGs
annotations.df <- annotations.df %>%
  filter(fasta %in% mags.df$MAG_ID)

# subset annotations to genes of interest
annotations.df <- annotations.df %>%
  filter(apply(annotations.df, 1, function(row) any(row %in% genes.df$Gene_ID)))

# subset columns 
annotations.df <- annotations.df %>% 
  select(fasta, ko_id, pfam_hits)

# edit pfam ids to remove trailing decimals
annotations.df <- annotations.df %>%
  mutate(pfam_hits = gsub("\\.(\\d+)(?=;|$)", "", pfam_hits, perl = TRUE))

# remove duplicated rows (some genomes have multiple copies of a gene)
annotations.df <- annotations.df[!duplicated(annotations.df), ]

# join dataframes 
annotations_j.df <- annotations.df %>%
  full_join(genes.df, by = c("ko_id" = "Gene_ID")) ## note: NAs are genes not found in any genomes

# group by mag to combine multiple genes
annotations_g.df <- annotations_j.df %>%
  group_by(Pathway, Total_genes, fasta) %>%
  summarise(ko_id = paste(unique(ko_id), collapse = "; "),
            Gene_name = paste(unique(Gene_name), collapse = "; ")) %>%
  ungroup()                 

# get number of genes found for each pathway in each MAG
annotations_g.df <- annotations_g.df %>%
  mutate(Genes_found = 1 + str_count(ko_id, ";")) 

# get proportion of pathway complete (cap at 100 due to redundant genes)
annotations_g.df <- annotations_g.df %>%
  mutate(Pathway_completion = pmin(Genes_found / Total_genes * 100, 100))

# edit column name and rearrange
annotations_g.df <- annotations_g.df %>%
  rename("MAG_id" = "fasta") %>%
  select(MAG_id, Pathway, Total_genes, Genes_found, ko_id, Gene_name, Pathway_completion)

## add results from dram distilled - better for pathway completion for carbon 
# and denitrification pathways due to redundant genes

# select cols
distilled_carbon.df <- distilled.df %>%
  select("genome", 
         "Citrate cycle (TCA cycle, Krebs cycle)", 
         "Reductive pentose phosphate cycle (Calvin cycle)")

distilled_denit.df <- distilled.df %>%
  select("genome", 
         "Nitrogen metabolism: nitrate => nitrite", 
         "Nitrogen metabolism: nitrite => nitric oxide",
         "Nitrogen metabolism: nitric oxide => nitrous oxide",
         "Nitrogen metabolism: nitrous oxide => nitrogen")

# sum the steps for dentifification (DRAM outputs T/F for each step; carbon already proportion of pathway)

# cols to sum
cols_to_sum <- colnames(distilled_denit.df)[-1]

# change from true /false character vector to numeric
distilled_denit.df[cols_to_sum] <- lapply(distilled_denit.df[cols_to_sum], function(x) {
  as.numeric(tolower(trimws(x)) == "true")
})

# sum rows
distilled_denit.df$total <- rowSums(distilled_denit.df[cols_to_sum], na.rm = TRUE)

# add pathway and pathway completion columns
distilled_denit.df$Pathway_completion <- distilled_denit.df$total / 4 * 100
distilled_denit.df$Pathway <- rep("Denitrification", nrow(distilled_denit.df))

# select columns
distilled_denit.df <- distilled_denit.df %>% 
  select(genome, Pathway_completion, Pathway) 

# rename cols
colnames(distilled_carbon.df) <- c("MAG_id", "TCA cycle", "Calvin cycle")
colnames(distilled_denit.df)[1] <- "MAG_id"

# change to long format
distilled_carbon.df <- distilled_carbon.df %>%
  pivot_longer(cols = c("TCA cycle", "Calvin cycle"), 
               names_to = "Pathway", 
               values_to = "Pathway_completion")

# change complete to %
distilled_carbon.df$Pathway_completion <- distilled_carbon.df$Pathway_completion * 100

# replace with new pathway completion values
# carbon
annotations_g.df <- annotations_g.df %>%
  left_join(distilled_carbon.df, 
            by = c("MAG_id", "Pathway"),
            suffix = c("", ".new")) %>%
  mutate(Pathway_completion = coalesce(Pathway_completion.new, Pathway_completion)) %>%
  select(-Pathway_completion.new)

# denitrification
annotations_g.df <- annotations_g.df %>%
  left_join(distilled_denit.df, 
            by = c("MAG_id", "Pathway"),
            suffix = c("", ".new")) %>%
  mutate(Pathway_completion = coalesce(Pathway_completion.new, Pathway_completion)) %>%
  select(-Pathway_completion.new)


## save output ----

write.table(annotations_g.df, 
            "Data/Metagenomics_full/Formatted_outputs/mags_genes_of_interest.tsv", 
            sep = "\t", row.names = F)

# import to script #15 to plot with tree

## combine with taxonomy ----

# import gene table
# note: when mag id = NA indicates genes not found in any genomes
genes.df <- read.table(
  "Data/Metagenomics_full/Formatted/mags_genes_of_interest.tsv",
  header = T, strip.white = T, sep = "\t"
)

# combine with taxonomy
genes.df <- left_join(genes.df, tax.df, join_by("MAG_id" == "MAG_ID"))

# combine mags with abundance to see which samples they are present in  
# get rownames and transpose
rownames(mags.df) <- mags.df$MAG_ID
mags.df <- t(mags.df[,-1]) %>% as.data.frame()

# get sample id
mags.df$Sample_ID <- rownames(mags.df)

# join
mags.df <- left_join(meta.df, mags.df)

# check taxa for different pathways
# photosystems
photo.df <- genes.df %>%
  filter(Pathway %in% c("Photosystem I (main subunits)", 
                        "Photosystem II (main subunits)"))

photo.df$Phylum %>% unique()

# calvin cycle
calvin.df <- genes.df %>%
  filter(Pathway == "Calvin cycle")

# carotenoids
beta.df <- genes.df %>%
  filter(Pathway == "Carotenoid (Beta-carotene pathway)")

beta_k.df <- genes.df %>%
  filter(str_detect(ko_id, "K06443"))

beta_k.df$Phylum %>% unique()

acid.df <- genes.df %>%
  filter(Pathway == "Carotenoid (Abscisic acid pathway)")

# TCA cycle
tca.df <- genes.df %>%
  filter(Pathway == "TCA cycle") 

# key genes (not in rTCA cycle)
tca_k.df <- tca.df %>%
  filter(grepl("K01647|K05942|K01659", ko_id)) %>% # key step 
  filter(Genes_found > 9)
 # filter(grepl("K01647|K05942|K01659|K00164|K00658|K00234|K00235|K00236|K00237|K25801|K00244|K00245|K00246|K00247", ko_id)) # can have alternative genes

# nitrate reduction
nr.df <- genes.df %>%
  filter(Pathway == "Assimilatory nitrate reduction")

nr.df$Phylum %>% unique()
nr.df$Family %>% unique()

# nitrogen fixation
nif.df <- genes.df %>%
  filter(Pathway == "Nitrogen fixation")

# glutamine synthesis
gmine.df <- genes.df %>%
  filter(Pathway == "Glutamine synthesis")

# glutamate synthesis
gmate.df <- genes.df %>%
  filter(Pathway == "Glutamate synthesis")

# denitrification
dn.df <- genes.df %>%
  filter(Pathway == "Denitrification")

# no synthesis
no.df <- genes.df %>%
  filter(Gene_name == "nos")

# arginine synthesis
ar.df <- genes.df %>%
  filter(Pathway == "Arginine synthesis") %>%
  filter(Genes_found == 2)

# gaba 
gaba.df <- genes.df %>%
  filter(Pathway == "GABA synthesis")

# type III SS
t3.df <- gaba.df <- genes.df %>%
  filter(Pathway == "type III SS")

# type II SS
t2.df <- genes.df %>%
  filter(Pathway == "type II SS") %>%
  filter(Pathway_completion > 50)

# sec-SRP
sec.df <- genes.df %>%
  filter(Pathway == "Sec-SRP") %>%
  filter(Pathway_completion > 50)

# flagella
flg.df <- genes.df %>%
  filter(Pathway == "Flagellar")

# MAGs without flagellar
noflg.df <- tax.df %>%
  filter(!MAG_ID %in% flg.df$MAG_id)

# lps
lps.df <- genes.df %>%
  filter(Pathway %in% c("Lipopolysaccharide biosynthesis (lipid A Raetz pathway)", 
                        "Lipopolysaccharide biosynthesis (lipid A Modification pathway)")) %>%
  filter(Pathway_completion > 50)

# modification pathway 
mod_lps.df <- lps.df %>%
  filter(grepl("lpxE|eptA|kdoH1|kdoH2|lpxF|lpxR", Gene_name))

# agmo
agmo.df <- genes.df %>%
  filter(Gene_name == "AGMO")


## check blast ----

# keep only strong matches
cyclo <- cyclo %>%
  filter(e.value < 1e-10, bit_score > 100)

tbp <- tbp %>%
  filter(e.value < 1e-10, bit_score > 100)

# remove node id from scaffold
cyclo$scaffold_id <- gsub("_NODE.*", "", cyclo$scaffold_id)
tbp$scaffold_id <- gsub("_NODE.*", "", tbp$scaffold_id)

# filter tbp results to the genes from pseudoalteromonas sp. PS5 
# blast search multiple genomes, this one is from Alker et al. 2023
tbp_f <- tbp %>%
  filter(grepl("KR011923", query_sequence_id)) %>% 
  select(query_sequence_id, scaffold_id)

# group by mag to combine multiple genes
tbp_f <- tbp_f %>%
  group_by(scaffold_id) %>%
  summarise(query_sequence_id = paste(unique(query_sequence_id), collapse = "; ")) %>%
  ungroup()

cyclo_f <- cyclo %>%
  group_by(scaffold_id) %>%
  summarise(query_sequence_id = paste(unique(query_sequence_id), collapse = "; ")) %>%
  ungroup()


# get number of genes in the biosynthesis pathway
tbp_f <- tbp_f %>%
  mutate(gene_number = 1 + str_count(query_sequence_id, ";")) 

cyclo_f <- cyclo_f %>%
  mutate(gene_number = 1 + str_count(query_sequence_id, ";")) 

# get proportion of pathway complete 
# 6 genes total for tbp
tbp_f$tbp <- (tbp_f$gene_number / 6) * 100

# 11 genes total for cyclo
cyclo_f$cyclo <- (cyclo_f$gene_number / 11) * 100

# remove MAGs without key cyclisation gene
cyclo_f <- cyclo_f %>%
  filter(grepl("PRUB680", query_sequence_id)) 

# join with tax
cyclo_f <- left_join(tax.df, cyclo_f, join_by(MAG_ID == scaffold_id))
tbp_f <- left_join(tax.df, tbp_f, join_by(MAG_ID == scaffold_id))

# check which mags have genes
cyc.df <- cyclo_f %>%
  filter(cyclo > 0)

tbp.df <- tbp_f %>%
  filter(tbp > 0)


## create table of genes in MAGs ----

# reload annotations df
annotations.df <- read.table(
  "Data/Metagenomics_full/Raw/dram/genomes/annotations_sparse.tsv",
  header = T, strip.white = T, sep = "\t", fill = TRUE
)

# reload genes of interest
genes.df <- read.csv(
  "Data/Metagenomics_full/Formatted/genes_of_interest.csv",
  header = T, strip.white = T
)

# subset annotations to genes of interest
annotations.df <- annotations.df %>%
  filter(apply(annotations.df, 1, function(row) any(row %in% genes.df$Gene_ID)))

# subset columns 
annotations.df <- annotations.df %>% 
  select(fasta, ko_id)

# join dataframes 
annotations_j.df <- annotations.df %>%
  full_join(genes.df, by = c("ko_id" = "Gene_ID")) %>%
  select(!Pathway)

# combine KO with gene name and edit columns
annotations_j.df <- annotations_j.df %>%
  mutate(Gene_name = paste0(Gene_name, " (", ko_id, ")")) %>%
  select(!ko_id) %>%
  rename(MAG_id = fasta)

# add blast results
# select cols
cyclo_f2 <- cyclo %>%
  filter(query_sequence_id == "PRUB680") %>%
  select(scaffold_id, query_sequence_id)

tbp_f2 <- tbp %>%
  select(scaffold_id, query_sequence_id)

# subset to pseudoalteromonas pS5 genome
tbp_f2 <- tbp_f2 %>%
  filter(grepl("KR011923", query_sequence_id))

# rename genes to gene names
tbp_f2$query_sequence_id <- gsub("lcl|KR011923.1_gene_1", "bmp10", tbp_f2$query_sequence_id, fixed = TRUE)
tbp_f2$query_sequence_id <- gsub("lcl|KR011923.1_gene_2", "bmp9", tbp_f2$query_sequence_id, fixed = TRUE)
tbp_f2$query_sequence_id <- gsub("lcl|KR011923.1_gene_4", "bmp2", tbp_f2$query_sequence_id, fixed = TRUE)
tbp_f2$query_sequence_id <- gsub("lcl|KR011923.1_gene_5", "bmp3", tbp_f2$query_sequence_id, fixed = TRUE)
tbp_f2$query_sequence_id <- gsub("lcl|KR011923.1_gene_6", "bmp4", tbp_f2$query_sequence_id, fixed = TRUE)

# edit colnames
colnames(cyclo_f2) <- c("MAG_id", "Gene_name")
colnames(tbp_f2) <- c("MAG_id", "Gene_name")

# add to annotations
annotations_j.df <- rbind(annotations_j.df, cyclo_f2)
annotations_j.df <- rbind(annotations_j.df, tbp_f2)

# make wide format and sum copy numbers for each gene
annotations_w.df <- annotations_j.df %>%
  count(MAG_id, Gene_name) %>%                      
  pivot_wider(names_from = Gene_name,           
              values_from = n,    
              values_fill = 0) 

# add taxonomy
annotations_w.df <- left_join(tax.df, annotations_w.df, join_by("MAG_ID" == "MAG_id"))

# save
write.table(annotations_w.df, "mags_with_genes_of_interest_supptable.tsv", 
            sep = "\t", row.names = F, col.names = T)



















