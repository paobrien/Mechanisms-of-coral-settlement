## plot metagenome phylogeny with genes of interest

setwd("~/Documents/R/Microbial_inducers")

library(tidyverse)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(scales)
library(viridis)
library(ggnewscale)


## import and format data ---- 

## formatted MAG data
mag.df <- read.table(
  "Data/Metagenomics_full/Formatted/dereplicated_mags_summary_data_filtered.tsv",
  sep = '\t', header = T, strip.white = T
)

# remove unmapped data
mag.df <- mag.df %>%
  filter(MAG_ID != "unmapped")

# reorder cols
mag.df <- mag.df %>%
  select(MAG_ID, starts_with("SE"), everything())

# change NAs to 0
mag.df[is.na(mag.df)] <- 0

## metadata
meta.df <- read.csv(
  "Data/Metagenomics_full/Formatted/sample_metadata_filtered.csv",
  header = T, strip.white = T
)

## GTDB tree
tree <- ape::read.tree(
  "Data/Metagenomics_full/Raw/gtdb/derep_bins.bac120.classify.tree"
)

# Subset tree to our bins
tree <- keep.tip(tree, mag.df$MAG_ID)
plot(tree, show.tip.label = F) # check tree

## import gene table
genes.df <- read.table(
  "Data/Metagenomics_full/Formatted/mags_genes_of_interest.tsv",
  header = T, strip.white = T, sep = "\t"
)

## import blast results
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

## plot tree ----

# edit tax labels
mag.df$Phylum <- gsub("p__", "", mag.df$Phylum)
mag.df$Class <- gsub("c__", "", mag.df$Class)
mag.df$Family <- gsub("f__$", "Unclassified", mag.df$Family)
mag.df$Family <- gsub("f__", "", mag.df$Family)
mag.df$Phylum <- mag.df$Phylum %>% as.factor()
mag.df$Class <- mag.df$Class %>% as.factor()
mag.df$Family <- mag.df$Family %>% as.factor()

# relevel factors 
mag.df <- mag.df[match(tree$tip.label, mag.df$MAG_ID),] # reorder rows to phylogeny
phylum <- mag.df$Phylum %>% unique

mag.df$Phylum <- factor(mag.df$Phylum, levels = phylum)
levels(mag.df$Phylum)

# get tree colours (discrete palette from chat)
# 20 colours
phy_col <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
  "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69"
)
show_col(phy_col)

# plot tree
p <- ggtree(tree, size = 0.25, layout = "fan", 
            open.angle = 20, ladderize = F)  %<+% mag.df +
  geom_tippoint(aes(fill = Phylum), 
                size = 2, 
                shape = 21,
                alpha = 1) +
  geom_tiplab(size = 0.000001, # (change to 1 if for visible tip labs - supp figure)
              align=TRUE, 
              linetype='dashed', 
              linesize=0.075, 
              offset = 0.35, 
              aes(label = Family))  +
  scale_fill_manual(values = phy_col) +
  guides(fill = guide_legend(ncol=1, reverse = T)) + 
  theme(legend.position="left",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 13))
p 

## add gene presence ----

## for blast results

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

# join with metadata
mag.df <- left_join(mag.df, cyclo_f, join_by(MAG_ID == scaffold_id))
mag.df <- left_join(mag.df, tbp_f, join_by(MAG_ID == scaffold_id))

# change NAs to 0
mag.df <- mag.df %>%
  mutate(across(c(cyclo, gene_number.x, tbp, gene_number.y), ~replace_na(., 0)))

## for KO results

# make wide format
genes_w.df <- genes.df %>%
  select(MAG_id, Pathway, Pathway_completion) %>%
  pivot_wider(names_from = c(Pathway), values_from = c(Pathway_completion))

# change NAs to 0
genes_w.df <- genes_w.df %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))

# remove NA mag id (genes that weren't in mags)
genes_w.df <- genes_w.df %>%
  filter(!is.na(MAG_id))

# add tbp and cyclo blast results
genes_w.df <- left_join(genes_w.df, select(mag.df, "MAG_ID", "cyclo", "tbp"), 
                        join_by("MAG_id" == "MAG_ID"))

# change back to long format 
genes_l.df <- genes_w.df %>%
  pivot_longer(cols = -c(MAG_id),
               names_to = "Pathway",
               values_to = "Pathway_completion")

# subset to inducing genes
genes_i_l.df <- genes_l.df %>%
  filter(Pathway %in% c("Assimilatory nitrate reduction",
                        "Carotenoid (Beta-carotene pathway)",
                        "Glutamate synthesis",
                        "GABA synthesis",
                        "type III SS",
                        "cyclo",
                        "tbp"))

# subset to carbon and nitrogen metabolism genes
genes_cn_l.df <- genes_l.df %>%
  filter(Pathway %in% c("Photosystem I (main subunits)",
                        "Photosystem II (main subunits)",
                        "Calvin cycle",
                        "TCA cycle",
                        "Assimilatory nitrate reduction",
                        "Denitrification",
                        "Nitrogen fixation",
                        "Nitric oxide synthesis"))

# subset to amino acids, neurotransmiter, SS, LPS, flagellar genes
genes_other_l.df <- genes_l.df %>%
  filter(Pathway %in% c("Glutamine synthesis",
                        "Glutamate synthesis",
                        "Arginine synthesis",
                        "GABA synthesis",
                        "Cycloprodigiosin biosynthesis",
                        "Sec-SRP",
                        "type II SS",
                        "type III SS",
                        "Flagellar",
                        "Lipopolysaccharide biosynthesis (lipid A Raetz pathway)"))

## add to tree plot

# inducing genes
# get custom x-axis labels
custom_labels <- c(
  "Assimilatory nitrate reduction" = "Nitrate reduction",
  "Carotenoid (Beta-carotene pathway)" = "Carotenoid biosynthesis",
  "Glutamate synthesis" = "Glutamate biosynthesis",
  "GABA synthesis" = "GABA biosynthesis",
  "type III SS" = "Type III SS",
  "cyclo" = "Cycloprodigiosin biosynthesis",
  "tbp" = "TBP biosynthesis"
)

# update factor labels
genes_i_l.df$Pathway <- factor(genes_i_l.df$Pathway, 
                               levels = names(custom_labels),
                               labels = custom_labels)

# plot
p + new_scale_fill() +
  geom_fruit(data = genes_i_l.df,
             geom = geom_tile, 
             mapping = aes(y=MAG_id, 
                           x=Pathway,
                           fill=Pathway_completion),
             offset = 0.02, # change to 0.01 for tiplabs
             pwidth = 0.11,
             axis.params = list(
               axis = "x",             
               text.size = 2.1,        
               text.angle = 90,        
               hjust = 1
             )) +
  scale_fill_viridis(option = "viridis",  
                     name = "Pathway completion (%)")

# save
ggsave("mag_tree_circular_with_inducers.svg", 
       device = "svg", 
       width = 34, 
       height = 34, 
       scale = 1, 
       units = "cm", 
       dpi = 300)

# carbon and nitrogen metabolism

# get custom x-axis labels
custom_labels <- c(
  "Photosystem I (main subunits)" = "Photosystem I",
  "Photosystem II (main subunits)" = "Photosystem II",
  "Calvin cycle" = "Calvin cycle",
  "TCA cycle" = "TCA cycle",
  "Assimilatory nitrate reduction" = "Nitrate reduction",
  "Denitrification" = "Denitrification",
  "Nitrogen fixation" = "Nitrogen fixation",
  "Nitric oxide synthesis" = "Nitric oxide synthesis (nos)"
)

# update factor labels
genes_cn_l.df$Pathway <- factor(genes_cn_l.df$Pathway, 
                               levels = names(custom_labels),
                               labels = custom_labels)

# redo tree for formatting
p <- ggtree(tree, size = 0.25, layout = "fan", 
            open.angle = 20, ladderize = F)  %<+% mag.df +
  geom_tippoint(aes(fill = Phylum), 
                size = 2, 
                shape = 21,
                alpha = 1) +
  geom_tiplab(size = 1, 
              align=TRUE, 
              linetype='dashed', 
              linesize=0.075, 
              offset = 0.4, 
              aes(label = Family))  +
  scale_fill_manual(values = phy_col) +
  guides(fill = guide_legend(ncol=1, reverse = T)) + 
  theme(legend.position="left",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 13))
p 

p + new_scale_fill() +
  geom_fruit(data = genes_cn_l.df,
             geom = geom_tile, 
             mapping = aes(y=MAG_id, 
                           x=Pathway,
                           fill=Pathway_completion),
             offset = 0.01, # change to 0.01 for tiplabs,
             pwidth = 0.125,
             axis.params = list(
               axis = "x",             
               text.size = 2.0,        
               text.angle = 90,        
               hjust = 1
             )) +
  scale_fill_viridis(option = "viridis",  
                     name = "Pathway completion (%)")

# save
ggsave("mag_tree_circular_carbon_nitrogen.svg", 
       device = "svg", 
       width = 34, 
       height = 34, 
       scale = 1, 
       units = "cm", 
       dpi = 300)

# amino acids, neurotransmitters, SS, flagellar, LPS

# get custom x-axis labels
custom_labels <- c(
  "Glutamine synthesis" = "Glutamine biosynthesis",
  "Glutamate synthesis" = "Glutamate biosynthesis",
  "Arginine synthesis" = "Arginine biosynthesis",
  "GABA synthesis" = "GABA biosynthesis",
  "Cycloprodigiosin biosynthesis" = "Alkylglycerol monooxygenase",
  "Sec-SRP" = "Sec-SRP pathway",
  "type II SS" = "type II SS",
  "type III SS" = "type III SS",
  "Flagellar" = "Flagella",
  "Lipopolysaccharide biosynthesis (lipid A Raetz pathway)" = "Lipopolysaccharide biosynthesis"
)

# update factor labels
genes_other_l.df$Pathway <- factor(genes_other_l.df$Pathway, 
                                levels = names(custom_labels),
                                labels = custom_labels)

# redo tree for formatting
p <- ggtree(tree, size = 0.25, layout = "fan", 
            open.angle = 20, ladderize = F)  %<+% mag.df +
  geom_tippoint(aes(fill = Phylum), 
                size = 2, 
                shape = 21,
                alpha = 1) +
  geom_tiplab(size = 1, 
              align=TRUE, 
              linetype='dashed', 
              linesize=0.075, 
              offset = 0.5, 
              aes(label = Family))  +
  scale_fill_manual(values = phy_col) +
  guides(fill = guide_legend(ncol=1, reverse = T)) + 
  theme(legend.position="left",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 13))
p 


p + new_scale_fill() +
  geom_fruit(data = genes_other_l.df,
             geom = geom_tile, 
             mapping = aes(y=MAG_id, 
                           x=Pathway,
                           fill=Pathway_completion),
             offset = 0.015,
             pwidth = 0.16,
             axis.params = list(
               axis = "x",             
               text.size = 2,        
               text.angle = 90,        
               hjust = 1
             )) +
  scale_fill_viridis(option = "viridis",  
                     name = "Pathway completion (%)")

# save
ggsave("mag_tree_circular_other_genes.svg", 
       device = "svg", 
       width = 35, 
       height = 35, 
       scale = 1, 
       units = "cm", 
       dpi = 300)




