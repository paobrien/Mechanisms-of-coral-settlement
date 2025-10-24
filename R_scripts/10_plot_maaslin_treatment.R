# format and plot maaslin output to check for genes that are differentially abundant among treatments

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)
library(KEGGREST)
library(ComplexHeatmap)
library(viridis)
library(scales)

## import and format data ----

# import significant results
maas.df <- read.table(
  "/Users/uqpobri2/Desktop/UQ_postdoc/Bioinformatics/MG_biofilm/10_maaslin/genes/treatment_no_ctrls/all_ko_lm_maaslin/significant_results.tsv", 
  header = T, 
  sep = "\t"
)

# ko frequency table
ko.df <- read.table(
  "Data/Metagenomics_full/Raw/dram/genes/ko_frequency_table.tsv", 
  header = T, sep = "\t", strip.white = T
)

# import metadata
meta.df <- read.csv(
  "Data/Metagenomics_full/Formatted/sample_metadata_filtered.csv",
  header = T, strip.white = T
)

# filter blanks and controls
meta.df <- meta.df %>%
  filter(Treatment != "negative") %>%
  filter(Treatment != "control")

# filter ko df
ko.df <- ko.df %>%
  dplyr::select(ko_id, meta.df$Sample_ID)

## format with KEGG rest ----

# remove uneeded columns
maas.df <- maas.df %>%
    dplyr::select(!c(metadata, N, N.not.0))

# get ko pathways
ko_pathways.l <- keggLink("ko", "pathway")

# create df
ko_pathways.df <- utils::stack(ko_pathways.l) %>% 
  dplyr::filter(!grepl("path:ko", ind)) %>%
  dplyr::mutate(values = gsub("ko:", "",values),
                ind = gsub("path:", "", ind)) %>%
  dplyr::rename("KO" = values, "Pathway" = ind)

# get pathway definitions
ko_pathway_defintions.l <- keggList("pathway") 

# create df
ko_pathway_defintions.df <- stack(ko_pathway_defintions.l)

# rename columns
names(ko_pathway_defintions.df) <- c("Pathway_description","Pathway")

# join with ko's
ko_pathways.df <- left_join(ko_pathways.df, ko_pathway_defintions.df)

# join with maaslin
ko_pathways.df <- left_join(maas.df, ko_pathways.df, by = c("feature" = "KO"))

# replace NA wth Unknown
ko_pathways.df$Pathway_description <- ifelse(
  is.na(ko_pathways.df$Pathway_description), "Unknown", ko_pathways.df$Pathway_description)

# subset to KOs more abundant in L_2M treatment (top 50)
ko_pos.df <- maas.df %>%
  filter(value == "D_2M" & coef < 0) %>%
  slice_head(n = 50)

ko_pos.df <- ko.df %>%
  filter(ko_id %in% ko_pos.df$feature)

# subset to KOs more abundant in D_2M treatment (top 50)
ko_neg.df <- maas.df %>%
  filter(value == "D_2M" & coef > 0) %>%
  slice_head(n = 50)

ko_neg.df <- ko.df %>%
  filter(ko_id %in% ko_neg.df$feature)

# get pathway label
path_lab.df <- ko_pathways.df %>%
  group_by(feature) %>%
  summarize(Pathway_description = paste(unique(Pathway_description), collapse = "/")) %>%
  ungroup() %>%
  mutate(label = paste0(feature, ": ", Pathway_description))

# join with with df
ko_pos.df <- left_join(ko_pos.df, path_lab.df, join_by(ko_id == feature))
ko_neg.df <- left_join(ko_neg.df, path_lab.df, join_by(ko_id == feature))

# edit long labels for plotting
# pos
ko_pos.df$label <- gsub("/Metabolic pathways/Biosynthesis of secondary metabolites", "", ko_pos.df$label)
ko_pos.df$label <- gsub("/Metabolic pathways/Microbial metabolism in diverse environments", "", ko_pos.df$label)
ko_pos.df$label <- gsub("/Metabolic pathways", "", ko_pos.df$label)
ko_pos.df$label <- gsub("/Biosynthesis of secondary metabolites", "", ko_pos.df$label)
ko_pos.df$label <- gsub("Microbial metabolism in diverse environments/", "", ko_pos.df$label)
ko_pos.df$label <- gsub("/Amyotrophic lateral sclerosis/Huntington disease/Pathways of neurodegeneration - multiple diseases", "", ko_pos.df$label)
ko_pos.df$label <- gsub("/Aflatoxin biosynthesis/Pyruvate metabolism/Propanoate metabolism/Fatty acid metabolism/AMPK signaling pathway/Insulin signaling pathway/Glucagon signaling pathway/Alcoholic liver disease", "", ko_pos.df$label)
ko_pos.df$label <- gsub("- Pseudomonas aeruginosa", "", ko_pos.df$label)
ko_pos.df$label <- gsub("/Biosynthesis of cofactors", "", ko_pos.df$label)
ko_pos.df$label <- gsub("/Chloroalkane and chloroalkene degradation/Naphthalene degradation/Butanoate metabolism/Degradation of aromatic compounds/Amoebiasis", "/Butanoate metabolism", ko_pos.df$label)

# neg
ko_neg.df$label <- gsub("/Metabolic pathways/Biosynthesis of secondary metabolites/Microbial metabolism in diverse environments", "", ko_neg.df$label)
ko_neg.df$label <- gsub("Metabolic pathways/Biosynthesis of secondary metabolites/", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Biosynthesis of secondary metabolites", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Metabolic pathways", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Microbial metabolism in diverse environments", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Carbon metabolism", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/2-Oxocarboxylic acid metabolism/Biosynthesis of amino acids", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/C5-Branched dibasic acid metabolism", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Biosynthesis of amino acids/HIF-1 signaling pathway/Alzheimer disease/Pathogenic Escherichia coli infection/Salmonella infection/Diabetic cardiomyopathy", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Lysine degradation/Tyrosine metabolism/Butanoate metabolism/Nicotinate and nicotinamide metabolism", "", ko_neg.df$label)
ko_neg.df$label <- gsub("beta-Alanine metabolism/Inositol phosphate metabolism/Propanoate metabolism", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/2-Oxocarboxylic acid metabolism", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Lipoic acid metabolism", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Propanoate metabolism", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Nitrogen metabolism", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Glycine, serine and threonine metabolism/Valine, leucine and isoleucine degradation/Lysine degradation/Tryptophan metabolism/Pyruvate metabolism/Glyoxylate and dicarboxylate metabolism/One carbon pool by folate/Biosynthesis of cofactors", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Thyroid hormone synthesis/Diabetic cardiomyopathy", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Fatty acid metabolism", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Biosynthesis of cofactors", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Other carbon fixation pathways", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Lysine degradation/Benzoate degradation/Tryptophan metabolism/beta-Alanine metabolism/Butanoate metabolism/Pinene, camphor and geraniol degradation/Caprolactam degradation", "", ko_neg.df$label)
ko_neg.df$label <- gsub(" - Pseudomonas aeruginosa", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Fructose and mannose metabolism/Galactose metabolism/Purine metabolism/Starch and sucrose metabolism/Amino sugar and nucleotide sugar metabolism/Streptomycin biosynthesis/Biosynthesis of nucleotide sugars", "", ko_neg.df$label)
ko_neg.df$label <- gsub(" - plant/Mineral absorption", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Mycolic acid biosynthesis", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/Glyoxylate and dicarboxylate metabolism", "", ko_neg.df$label)
ko_neg.df$label <- gsub("/$"  , "", ko_neg.df$label)

# get grouping label
ko_neg.df$group <- ko_neg.df$label
ko_neg.df$group <- gsub("^K\\d+:\\s*", "", ko_neg.df$group) 
ko_neg.df$group <- gsub("/Carbon fixation by Calvin cycle", "", ko_neg.df$group) 
ko_neg.df$group <- gsub("Glycolysis / Gluconeogenesis/Citrate cycle (TCA cycle)", "Citrate cycle (TCA cycle)", ko_neg.df$group, fixed = T) 
ko_neg.df$group <- gsub("/Alanine, aspartate and glutamate metabolism", "", ko_neg.df$group) 
ko_neg.df$group <- gsub("/Lysine biosynthesis", "", ko_neg.df$group)
ko_neg.df$group <- gsub("/Glyoxylate and dicarboxylate metabolism", "", ko_neg.df$group) 
ko_neg.df$group <- gsub("Fatty acid degradation/", "", ko_neg.df$group) 
ko_neg.df$group <- gsub("Fatty acid biosynthesis/", "", ko_neg.df$group) 
ko_neg.df$group <- gsub("/Pentose phosphate pathway", "", ko_neg.df$group) 
ko_neg.df$group <- gsub("Valine, leucine and isoleucine biosynthesis", "Valine, leucine and isoleucine degradation", ko_neg.df$group) 


## plot heatmap ----

# make matrix
rownames(ko_pos.df) <- ko_pos.df$label
ko_pos.mat <- ko_pos.df %>%
  dplyr::select(!c(ko_id, Pathway_description, label))

rownames(ko_neg.df) <- ko_neg.df$label
ko_neg.mat <- ko_neg.df %>%
  dplyr::select(!c(ko_id, Pathway_description, label, group))

# format metadata
rownames(meta.df) <- meta.df$Sample_ID

# format metadata labels
meta.df$Treatment <- factor(meta.df$Treatment,
                            levels = c("D_2M", "L_1M", "L_2M"),
                            labels = c("Dark 2M", "Light 1M", "Light 2M"))

# select variables
meta_ha <- meta.df %>% 
  dplyr::select(Treatment)

# get column annotations
ha <- HeatmapAnnotation(
  df = meta_ha,
  col = list(
    Treatment = c("Dark 2M" = "#30123BFF", 
                  "Light 1M" = "#1AE4B6FF", 
                  "Light 2M" = "#FABA39FF"))
)

# plot heatmap
# pos
heat_cols <- viridis(10, option = "turbo")
ht <- Heatmap(ko_pos.mat, 
              name = "nRPKM",
              col = heat_cols,
              top_annotation = ha, 
              show_column_names = FALSE,
              row_names_gp = gpar(fontsize = 12),
              row_names_max_width = unit(15, "cm"),
              row_split = ko_pos.df$Pathway_description,
              row_dend_reorder = TRUE,
              column_dend_reorder = TRUE,
              row_title = NULL,
              column_split = factor(as.character(meta.df$Treatment),
                                    levels = c("Light 2M", "Light 1M", "Dark 2M")),
              border = TRUE)
draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

# save
svg("gene_treatment_pos_heatmap.svg", width = 13, height = 10)
draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")
dev.off()


pdf("gene_treatment_pos_heatmap.pdf", width = 13, height = 10)

draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

dev.off()


# neg
ht <- Heatmap(ko_neg.mat, 
              name = "nRPKM",
              col = heat_cols,
              top_annotation = ha, 
              show_column_names = FALSE,
              row_names_gp = gpar(fontsize = 12),
              row_names_max_width = unit(15, "cm"),
              row_split = ko_neg.df$group,
              row_dend_reorder = TRUE,
              column_dend_reorder = TRUE,
              row_title = NULL,
              column_split = factor(as.character(meta.df$Treatment),
                                    levels = c("Dark 2M", "Light 2M", "Light 1M")),
              border = TRUE)
draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

# to save
svg("gene_treatment_neg_heatmap.svg", width = 13, height = 10)
draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")
dev.off()


pdf("gene_treatment_neg_heatmap.pdf", width = 13, height = 10)
draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

dev.off()

