## plot figures and summarise genes of interest

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(scales)
library(ggpubr)

## import data ----

# ko frequency table
ko.df <- read.table(
  "Data/Metagenomics_full/Raw/dram/genes/ko_frequency_table.tsv",
  sep = "\t", header = T, strip.white = T
)

# load metadata
meta.df <- read.csv(
  "Data/Metagenomics_full/Formatted/sample_metadata_filtered.csv",
  header = T, strip.white = T
)

# genes of interest
genes.df <- read.csv(
  "Data/Metagenomics_full/Formatted/genes_of_interest.csv",
  header = T, strip.white = T
)

# load maaslin results

# settlement model
# get file path
fp <- "/Users/uqpobri2/Desktop/UQ_postdoc/Bioinformatics/MG_biofilm/10_maaslin/genes/no_ctrl/" 

# get directory list
dirs <- list.files(fp, pattern = "_maaslin")

# import significant results from each folder
maas_set.l <- list()
for (i in seq_along(dirs)) {
  maas_set.l[[i]] <- read.table(paste0(fp, dirs[i], "/significant_results.tsv"),
                               header = T,
                               sep = "\t")
  names(maas_set.l)[i] <- gsub("_ko_lm_maaslin", "", dirs[i])
}

# treatment model
# import significant results
maas_trmt.df <- read.table(
  "/Users/uqpobri2/Desktop/UQ_postdoc/Bioinformatics/MG_biofilm/10_maaslin/genes/treatment_no_ctrls/all_ko_lm_maaslin/significant_results.tsv",
  header = T,
  sep = "\t"
  )


## format data ----

# remove blanks and controls
meta.df <- meta.df %>%
  filter(Treatment != "negative") %>%
  filter(Treatment != "control")

# format metadata labels
meta.df$Treatment <- factor(meta.df$Treatment,
                            levels = c("D_2M", "L_1M", "L_2M"),
                            labels = c("Dark 2M", "Light 1M", "Light 2M"))

meta.df$Coral <- factor(meta.df$Coral,
                        levels = c("P_sinensis", "D_favus", "E_aspera", "P_lobata"),
                        labels = c("P. sinensis", "D. favus", "E. aspera", "P. lobata"))

ko.df <- ko.df %>%
  select(ko_id, meta.df$Sample_ID)

## add maaslin results

# combine to one df
# first add coral column
maas_set.l <- lapply(names(maas_set.l), function(name) {
  df <- maas_set.l[[name]]
  df$Coral <- rep(name, nrow(df))
  return(df)
})
maas_trmt.df$Coral <- rep("Treatment", nrow(maas_trmt.df))

# combine
maas.df <- do.call(rbind, maas_set.l)
maas.df <- rbind(maas.df, maas_trmt.df)

# select column
maas.df <- maas.df %>%
  select(!c("metadata", "N", "N.not.0"))

# join with genes list
genes_j.df <- left_join(genes.df, maas.df, join_by("Gene_ID" == "feature"))

# group by KO
genes_g.df <- genes_j.df %>%
  group_by(Gene_ID, Gene_name, Pathway) %>%
  summarise(
    value = paste(value, collapse = "; "),
    coef = paste(coef, collapse = "; "),
    pval = paste(pval, collapse = "; "),
    qval = paste(qval, collapse = "; "),
    Coral = paste(Coral, collapse = "; ")
)
genes_g.df # majority of genes are signficant

## save table
write.table(genes_j.df, file = "genes_of_interest_results.tsv",
            sep = "\t", row.names = F, col.names = T)

## plot - treatment ----

# subset genes associated with treatment
# all treatment
genes_treatment.df <- genes.df %>%
  filter(Pathway %in% c("Photosystem II (main subunits)", 
                        "Photosystem I (main subunits)", 
                        "Calvin cycle", 
                        "TCA cycle", 
                        "Carotenoid (Beta-carotene pathway)", 
                        "Carotenoid (Abscisic acid pathway)",
                        "Assimilatory nitrate reduction", 
                        "Denitrification", 
                        "Nitrogen fixation",
                        "Glutamine synthesis", 
                        "Glutamate synthesis", 
                        "GABA synthesis", 
                        "Nitric oxide synthesis",
                        "Nitric oxide regulation",
                        "Arginine synthesis")) 

ko_treatment.df <- ko.df %>%
  filter(ko_id %in% genes_treatment.df$Gene_ID)

# create row label 
genes_treatment.df$Pathway <- gsub("\\(main subunits\\)", "", genes_treatment.df$Pathway)
genes_treatment.df <- genes_treatment.df %>%
  mutate(label = paste0(Pathway, " (", Gene_name, ")"))

ko_treatment.df <- left_join(
  ko_treatment.df, genes_treatment.df, join_by("ko_id" == "Gene_ID")
  )

ko_treatment.df$label <- gsub("(Beta-carotene pathway)", "biosynthesis", ko_treatment.df$label, fixed = T)

# make matrix
rownames(ko_treatment.df) <- ko_treatment.df$label
ko_treatment.mat <- ko_treatment.df %>%
  select(!c(ko_id, Gene_name, Pathway, label))

# transform matrix (hard to see patterns with different gene abundances)
ko_treatment_log10.mat <- log10(ko_treatment.mat + 1)

# format metadata
rownames(meta.df) <- meta.df$Sample_ID

## subset to plot different pathways

# nitrogen metabolism
ko_treatment_nit.df <- ko_treatment.df %>%
  filter(Pathway %in% c("Assimilatory nitrate reduction",
                        "Denitrification", 
                        "Nitrogen fixation"))

ko_treatment_nit_log10.mat <- ko_treatment_log10.mat %>%
  mutate(Gene = rownames(ko_treatment_log10.mat)) %>%
  filter(Gene %in% ko_treatment_nit.df$label) %>%
  select(!Gene)

# carbon metabolism 
ko_treatment_car.df <- ko_treatment.df %>%
  filter(Pathway %in% c("Calvin cycle", "TCA cycle"))

ko_treatment_car_log10.mat <- ko_treatment_log10.mat %>%
  mutate(Gene = rownames(ko_treatment_log10.mat)) %>%
  filter(Gene %in% ko_treatment_car.df$label) %>%
  select(!Gene)

# photosystems and carotenoids
ko_treatment_ph.df <- ko_treatment.df %>%
  filter(Pathway %in% c("Photosystem I ", 
                        "Photosystem II ", 
                        "Carotenoid (Beta-carotene pathway)"))

ko_treatment_ph_log10.mat <- ko_treatment_log10.mat %>%
  mutate(Gene = rownames(ko_treatment_log10.mat)) %>%
  filter(Gene %in% ko_treatment_ph.df$label) %>%
  select(!Gene)

# amino acids and NO
ko_treatment_aa.df <- ko_treatment.df %>%
  filter(Pathway %in% c("Assimilatory nitrate reduction", 
                        "Glutamate synthesis", 
                        "Glutamine synthesis", 
                        "Arginine synthesis",
                        "GABA synthesis",
                        "Nitric oxide synthesis",
                        "Nitric oxide regulation"))

ko_treatment_aa_log10.mat <- ko_treatment_log10.mat %>%
  mutate(Gene = rownames(ko_treatment_log10.mat)) %>%
  filter(Gene %in% ko_treatment_aa.df$label) %>%
  select(!Gene)

## plot heatmap

# select variables
meta_ha <- meta.df %>% 
  select(Treatment)

# get column annotations
ha <- HeatmapAnnotation(
  df = meta_ha,
  col = list(
    Treatment = c("Dark 2M" = "#30123BFF", 
                  "Light 1M" = "#1AE4B6FF", 
                  "Light 2M" = "#FABA39FF")),
  annotation_name_gp = gpar(fontsize = 11)  
)

# photosystems
heat_cols <- viridis(10, option = "turbo")
ht <- Heatmap(ko_treatment_ph_log10.mat, 
              name = "nRPKM (log)",
              col = heat_cols,
              top_annotation = ha,
              show_column_names = FALSE,
              row_names_gp = gpar(fontsize = 10),
              row_names_max_width = unit(15, "cm"),
              row_split = ko_treatment_ph.df$Pathway,
              row_dend_reorder = TRUE,
              column_dend_reorder = TRUE,
              row_title = NULL,
              column_split = factor(as.character(meta.df$Treatment),
                                    levels = c("Light 1M",
                                               "Light 2M", 
                                               "Dark 2M")),
              column_title_gp = gpar(fontsize = 12),
              border = TRUE)
draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

# save
pdf("photosystems_treatment.pdf", width = 9.7, height = 3.5)
svg("photosystems_treatment.svg", width = 9.7, height = 3.5)

draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

dev.off()

# nitrogen
ht <- Heatmap(ko_treatment_nit_log10.mat, 
              name = "nRPKM (log)",
              col = heat_cols,
              top_annotation = ha,
              show_column_names = FALSE,
              row_names_gp = gpar(fontsize = 10),
              row_names_max_width = unit(10, "cm"),
              row_split = ko_treatment_nit.df$Pathway,
              row_dend_reorder = TRUE,
              column_dend_reorder = TRUE,
              row_title = NULL,
              column_split = factor(meta.df$Treatment),
              column_title_gp = gpar(fontsize = 12),
              border = TRUE)
draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

pdf("nitrogen_metabolism_treatment_no_nr.pdf", width = 9, height = 3)
svg("nitrogen_metabolism_treatment_no_nr.svg", width = 9, height = 3)

draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

dev.off()


## plot single genes of interest with box plots 
aa.df <- ko_treatment_aa.df %>%
  pivot_longer(starts_with("SE"), names_to = "Sample_ID", values_to = "nRPKM")

# group nosip and nos into same pathway
aa.df$Pathway <- gsub("regulation", "synthesis", aa.df$Pathway)

# edit gene name
aa.df$Gene_name <- gsub("GLUD1_2, gdhA", "GLUD1_2", aa.df$Gene_name)

# attach metadata
aa.df <- right_join(meta.df, aa.df)

# change factor levels
aa.df$Pathway <- factor(aa.df$Pathway,
                         levels = c("Assimilatory nitrate reduction",
                                    "Glutamine synthesis", 
                                    "Glutamate synthesis",
                                    "GABA synthesis",
                                    "Nitric oxide synthesis", 
                                    "Arginine synthesis"))

# plot boxplot
p <- ggplot(aa.df, aes(x = Gene_name, y = nRPKM, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("Light 2M" = "#FABA39FF", 
                               "Light 1M" = "#1AE4B6FF", 
                               "Dark 2M" = "#30123BFF")) +
  scale_y_continuous(limits = c(0, NA)) +
  ylab("Gene abundance (nRPKM)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11), 
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 13),
        strip.text = element_text(size = 12))
p + facet_wrap("Pathway", scales = "free")

# save
ggsave("aa_treatment_boxplot.svg", 
       device = "svg", 
       width = 20, 
       height = 16, 
       units = "cm", 
       dpi = 300)

## summarise data - treatment ----

# function get mean and SE of for each group for a particular variable
# take a dataframe, grouping variables and response variable as input
# can take 1 or 2 groups
sum_stats <- function(df, group1, group2 = NULL, variable) {
  # format  columns (convert string to symbol)
  group1_col <- sym(group1)
  variable_col <- sym(variable)
  
  # if no group2 (NULL)
  if (is.null(group2)) {
    new_df <- df %>%
      group_by(!!group1_col) %>% 
      summarise(mean = mean(!!variable_col, na.rm = TRUE), 
                stdev = sd(!!variable_col, na.rm = TRUE), 
                sterr = sd(!!variable_col, na.rm = TRUE)/sqrt(n()),
                min = min(!!variable_col, na.rm = TRUE),
                max = max(!!variable_col, na.rm = TRUE))
  } else {
    group2_col <- sym(group2)
    # if group2   
    new_df <- df %>%
      group_by(!!group1_col, !!group2_col) %>% 
      summarise(mean = mean(!!variable_col, na.rm = TRUE), 
                stdev = sd(!!variable_col, na.rm = TRUE), 
                sterr = sd(!!variable_col, na.rm = TRUE)/sqrt(n()),
                min = min(!!variable_col, na.rm = TRUE),
                max = max(!!variable_col, na.rm = TRUE))
  }
  
  return(new_df)
}

# photosystems/carotenoid
# make data long format
ko_treatment_ph_l.df <- ko_treatment_ph.df %>%
  pivot_longer(cols = starts_with("SE"),
               names_to = "Sample_ID",
               values_to = "Abundance")

# add metadata
ko_treatment_ph_l.df <- right_join(meta.df, ko_treatment_ph_l.df)

# summarise
photo_sum <- sum_stats(ko_treatment_ph_l.df, 
                       group1 = "Pathway", 
                       group2 = "Treatment",
                       variable = "Abundance")
# round
photo_sum <- photo_sum %>%
  mutate(across(where(is.numeric), ~ round(., 3)))
photo_sum

# carbon metabolism
# make data long format
ko_treatment_car_l.df <- ko_treatment_car.df %>%
  pivot_longer(cols = starts_with("SE"),
               names_to = "Sample_ID",
               values_to = "Abundance")

# add metadata
ko_treatment_car_l.df <- right_join(meta.df, ko_treatment_car_l.df)

# summarise
carbon_sum <- sum_stats(ko_treatment_car_l.df, 
                        group1 = "Pathway", 
                        group2 = "Treatment",
                        variable = "Abundance")
carbon_sum

# nitrogen metabolism
# make data long format
ko_treatment_nit_l.df <- ko_treatment_nit.df %>%
  pivot_longer(cols = starts_with("SE"),
               names_to = "Sample_ID",
               values_to = "Abundance")

# add metadata
ko_treatment_nit_l.df <- right_join(meta.df, ko_treatment_nit_l.df)

# summarise
nitrogen_sum <- sum_stats(ko_treatment_nit_l.df, 
                          group1 = "Pathway", 
                          group2 = "Treatment",
                          variable = "Abundance")
nitrogen_sum

# amino acid metabolism
# summarise (already formatted above )
aa_sum <- sum_stats(aa.df, 
                    group1 = "Gene_name", 
                    group2 = "Treatment",
                    variable = "nRPKM")
aa_sum

# round
aa_sum <- aa_sum %>%
  mutate(across(where(is.numeric), ~ round(., 5)))


## plot heatmap - settlement ----

# subset to genes of interest
# all settlement
genes_settlement.df <- genes.df %>%
  filter(Pathway %in% c("Cycloprodigiosin biosynthesis",
                        "Sec-SRP",
                        "type II SS",
                        "type III SS",
                        "GABA synthesis",
                        "Flagellar",
                        "Lipopolysaccharide biosynthesis (lipid A Raetz pathway)")) 

ko_settlement.df <- ko.df %>%
  filter(ko_id %in% genes_settlement.df$Gene_ID)

# remove dark treatment
to_rm <- meta.df %>%
  filter(Treatment == "Dark 2M") %>%
  select(Sample_ID)

ko_settlement.df <- ko_settlement.df %>%
  select(!to_rm$Sample_ID)

meta_light <- meta.df %>%
  filter(!Sample_ID %in% to_rm$Sample_ID)

# create row label 
genes_settlement.df$Pathway <- gsub("\\(lipid A Raetz pathway\\)", "", genes_settlement.df$Pathway)
genes_settlement.df$Pathway <- gsub("Lipopolysaccharide biosynthesis ", "LPS biosynthesis", genes_settlement.df$Pathway)
genes_settlement.df <- genes_settlement.df %>%
  mutate(label = paste0(Pathway, " (", Gene_name, ")"))

ko_settlement.df <- left_join(
  ko_settlement.df, genes_settlement.df, join_by("ko_id" == "Gene_ID")
)

# make matrix
rownames(ko_settlement.df) <- ko_settlement.df$label
ko_settlement.mat <- ko_settlement.df %>%
  select(!c(ko_id, Gene_name, Pathway, label))

# type II + sec
ko_t2_sec.df <- ko_settlement.df %>%
  filter(Pathway %in% c("type II SS", 
                        "Sec-SRP"))

ko_t2_sec.mat <- ko_settlement.mat %>%
  mutate(Gene = rownames(ko_settlement.mat)) %>%
  filter(Gene %in% ko_t2_sec.df$label) %>%
  select(!Gene)

ko_t2_sec.mat <- log(ko_t2_sec.mat + 1)

# type III
ko_t3.df <- ko_settlement.df %>%
  filter(Pathway %in% c("type III SS"))

ko_t3.mat <- ko_settlement.mat %>%
  mutate(Gene = rownames(ko_settlement.mat)) %>%
  filter(Gene %in% ko_t3.df$label) %>%
  select(!Gene)

ko_t3.mat <- log(ko_t3.mat + 1)

## lps and flagella
ko_fl.df <- ko_settlement.df %>%
  filter(Pathway %in% c("Flagellar",
                        "LPS biosynthesis")) 

ko_fl.mat <- ko_settlement.mat %>%
  mutate(Gene = rownames(ko_settlement.mat)) %>%
  filter(Gene %in% ko_fl.df$label) %>%
  select(!Gene)

ko_fl.mat <- log(ko_fl.mat + 1)

## plot

# select variables
meta_ha <- meta_light %>% 
  select(Coral, Percent_settled)

# change to Settlement (%)
colnames(meta_ha) <- c("Coral", "Settlement (%)")

# get colours
# settlement
clrs <- viridis(8, option = "plasma")
show_col(clrs)
s_clrs <- colorRamp2(c(0, 50, 100), c("#0D0887FF", "#B93289FF", "#F0F921FF"))

# coral (to pick colours below)
c_clrs <- viridis(16, option = "mako")
show_col(c_clrs)
c_clrs

# heatmap
heat_cols <- viridis(10, option = "turbo")

# get column annotations
ha <- HeatmapAnnotation(
  df = meta_ha,
  col = list(
    Coral = c("P. sinensis" = "#60CEACFF",
              "D. favus" = "#3484A5FF",
              "E. aspera" = "#403872FF",
              "P. lobata" = "#0B0405FF"),
    `Settlement (%)` = s_clrs),
  annotation_legend_param = list(
    Coral = list(labels_gp = gpar(fontface = "italic", fontsize = 10)),
    `Settlement (%)` = list(labels_gp = gpar(fontface = "plain", fontsize = 10)))
)


# plot type 2 + sec
ht <- Heatmap(ko_t2_sec.mat, 
              name = "nRPKM (log)",
              col = heat_cols,
              top_annotation = ha, 
              show_column_names = FALSE,
              row_names_gp = gpar(fontsize = 10),
              row_names_max_width = unit(15, "cm"),
              row_split = factor(as.character(ko_t2_sec.df$Pathway)),
              clustering_distance_rows = "euclidean",
              row_dend_reorder = TRUE,
              column_dend_reorder = TRUE,
              row_title = NULL,
              column_split = factor(as.character(meta_light$Coral), 
                                    levels = c("E. aspera",
                                               "P. lobata",
                                               "D. favus",
                                               "P. sinensis")),
              clustering_distance_columns = "euclidean",
              border = TRUE,
              column_title_gp = gpar(fontface = "italic", fontsize = 12))
draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

# t2 + sec
pdf("type2_sec_heatmap.pdf", width = 11, height = 4.2)
svg("type2_sec_heatmap.svg", width = 11, height = 4.2)

draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

dev.off()

# plot type 3
ht <- Heatmap(ko_t3.mat, 
              name = "nRPKM (log)",
              col = heat_cols,
              top_annotation = ha, 
              show_column_names = FALSE,
              row_names_gp = gpar(fontsize = 10),
              row_names_max_width = unit(15, "cm"),
              row_split = factor(as.character(ko_t3.df$Pathway)),
              clustering_distance_rows = "euclidean",
              row_dend_reorder = TRUE,
              column_dend_reorder = TRUE,
              row_title = NULL,
              column_split = factor(as.character(meta_light$Coral), 
                                    levels = c("P. sinensis",
                                               "D. favus",
                                               "E. aspera",
                                               "P. lobata")),
              clustering_distance_columns = "euclidean",
              border = TRUE,
              column_title_gp = gpar(fontface = "italic", fontsize = 12))
draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

# save
pdf("type3_SS_heatmap.pdf", width = 11, height = 3.1)
svg("type3_SS_heatmap.svg", width = 11, height = 3.1)

draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

dev.off()


# plot lps + flagellar
ht <- Heatmap(ko_fl.mat, 
              name = "nRPKM (log)",
              col = heat_cols,
              top_annotation = ha, 
              show_column_names = FALSE,
              row_names_gp = gpar(fontsize = 10),
              row_names_max_width = unit(15, "cm"),
              row_split = factor(as.character(ko_fl.df$Pathway),
                                 levels = c("LPS biosynthesis", "Flagellar")),
              clustering_distance_rows = "euclidean",
              row_dend_reorder = TRUE,
              column_dend_reorder = TRUE,
              row_title = NULL,
              column_split = factor(as.character(meta_light$Coral), 
                                    levels = c("P. sinensis",
                                               "D. favus",
                                               "E. aspera",
                                               "P. lobata")),
              clustering_distance_columns = "euclidean",
              border = TRUE,
              column_title_gp = gpar(fontface = "italic", fontsize = 12))
draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

# save
pdf("lps_fl_heatmap.pdf", width = 11, height = 4.1)
svg("lps_fl_heatmap.svg", width = 11, height = 4.1)

draw(ht, 
     annotation_legend_side = "left", 
     heatmap_legend_side = "left")

dev.off()


## make regression plots for single genes of interest: GABA and AGMO

# subset genes
gaba.df <- ko_settlement.df %>%
  filter(Gene_name %in% c("gadAB"))

agmo.df <- ko_settlement.df %>%
  filter(Gene_name %in% c("AGMO"))

# change to long format
gaba.df <- gaba.df %>%
  pivot_longer(starts_with("SE"), names_to = "Sample_ID", values_to = "nRPKM")

agmo.df <- agmo.df %>%
  pivot_longer(starts_with("SE"), names_to = "Sample_ID", values_to = "nRPKM")

# attach metadata
gaba.df <- right_join(meta.df, gaba.df)
agmo.df <- right_join(meta.df, agmo.df)

# get facet labels
#f_labs <- c(
#  "D_favus" = "D. favus", 
#  "E_aspera" = "E. aspera",
#  "P_lobata" = "P. lobata",
#  "P_sinensis" = "P. sinensis"
#)

# plot gaba
p <- ggplot(gaba.df, aes(x = Percent_settled, y = nRPKM, colour = Coral)) +
  geom_point() + geom_smooth(method = "lm") +
  scale_colour_manual(values =  c("P. lobata" = "#0B0405FF", 
                                  "E. aspera" = "#403872FF",
                                  "D. favus" = "#3484A5FF",
                                  "P. sinensis" = "#60CEACFF")) +
  xlab("Settlement (%)") +
  ylab("gadAB abundance (nRPKM)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 11), 
       axis.title.x = element_text(size = 13),
       axis.text.y = element_text(size = 11), 
       axis.title.y = element_text(size = 13),
       legend.position = "none",
       strip.text = element_text(size = 12, face = "italic")) +
  facet_wrap("Coral")#, labeller = labeller(Coral = f_labs)) #+
#  stat_cor(aes(label=paste(after_stat(r.label),
#                            after_stat(p.label),
#                            sep ="~`;`~"),
#               colour="black"),
#           method="pearson",
#           size=4,
#           label.y=0.0235,
#           label.x = 0, digits = 3, p.accuracy = 0.001, label.sep = "; ")
p

# plot agmo
p2 <- ggplot(agmo.df, aes(x = Percent_settled, y = nRPKM, colour = Coral)) +
  geom_point() + geom_smooth(method = "lm") +
  scale_colour_manual(values =  c("P. lobata" = "#0B0405FF", 
                                  "E. aspera" = "#403872FF",
                                  "D. favus" = "#3484A5FF",
                                  "P. sinensis" = "#60CEACFF")) +
  xlab("Settlement (%)") +
  ylab("gadAB abundance (nRPKM)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 11), 
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 11), 
        axis.title.y = element_text(size = 13),
        legend.position = "none",
        strip.text = element_text(size = 12, face = "italic")) +
  facet_wrap("Coral")#, labeller = labeller(Coral = f_labs)) #+
#  stat_cor(aes(label=paste(after_stat(r.label),
#                            after_stat(p.label),
#                            sep ="~`;`~"),
#               colour="black"),
#           method="pearson",
#           size=4,
#           label.y=0.0235,
#           label.x = 0, digits = 3, p.accuracy = 0.001, label.sep = "; ")
p2

# run pearson correlations

# log transform using maaslin log function
LOG <- function(x) {
  y <- replace(x, x == 0, min(x[x>0]) / 2)
  return(log2(y))
}

# appy to df
ko_settlement_log.df <- ko_settlement.df %>%
  mutate(across(where(is.numeric), LOG))

# subset gaba and agmo
gaba_log.df <- ko_settlement_log.df %>%
  filter(Gene_name %in% c("gadAB"))

agmo_log.df <- ko_settlement_log.df %>%
  filter(Gene_name %in% c("AGMO"))

# change to long format
gaba_log.df <- gaba_log.df %>%
  pivot_longer(starts_with("SE"), names_to = "Sample_ID", values_to = "nRPKM")

agmo_log.df <- agmo_log.df %>%
  pivot_longer(starts_with("SE"), names_to = "Sample_ID", values_to = "nRPKM")

# attach metadata
gaba_log.df <- right_join(meta.df, gaba_log.df)
agmo_log.df <- right_join(meta.df, agmo_log.df)

# gaba
#psin
p_psin <- with(gaba_log.df %>% filter(., Coral == "P. sinensis"), cor.test(Percent_settled, nRPKM))$p.value
r_psin <- with(gaba_log.df %>% filter(., Coral == "P. sinensis"), cor.test(Percent_settled, nRPKM))$estimate

#dfav
p_dfav <- with(gaba_log.df %>% filter(., Coral == "D. favus"), cor.test(Percent_settled, nRPKM))$p.value
r_dfav <- with(gaba_log.df %>% filter(., Coral == "D. favus"), cor.test(Percent_settled, nRPKM))$estimate

#plob
p_plob <- with(gaba_log.df %>% filter(., Coral == "P. lobata"), cor.test(Percent_settled, nRPKM))$p.value
r_plob <- with(gaba_log.df %>% filter(., Coral == "P. lobata"), cor.test(Percent_settled, nRPKM))$estimate

#easp
p_easp <- with(gaba_log.df %>% filter(., Coral == "E. aspera"), cor.test(Percent_settled, nRPKM))$p.value
r_easp <- with(gaba_log.df %>% filter(., Coral == "E. aspera"), cor.test(Percent_settled, nRPKM))$estimate

# add annotations (adding stat cor manually due to formatting issues)
panel_labels <- data.frame(
  Coral = c("P. sinensis", "D. favus", "E. aspera", "P. lobata"),
  label = c("Coef = 0.142; FDR = 0.136\nR = 0.523; p < 0.001",
            "Coef = 0.060; FDR = 0.736\nR = 0.611; p < 0.001", 
            "Coef = -0.181; FDR = 0.059\nR = 0.111; p = 0.148", 
            "Coef = 0.201; FDR = 0.216\nR = 0.447; p = 0.003"),
  Percent_settled = 0,   
  nRPKM = 0.0248
)

panel_labels$Coral <- factor(panel_labels$Coral,
                             levels = levels(gaba.df$Coral))

p +
  geom_text(data = panel_labels,
            aes(x = Percent_settled, y = nRPKM, hjust = 0, label = label),
            inherit.aes = FALSE, size = 3.5)

# save
ggsave("gaba_regression_plot.svg", 
       device = "svg", 
       width = 20, 
       height = 14, 
       units = "cm", 
       dpi = 300)

#agmo
#psin
p_psin <- with(agmo_log.df %>% filter(., Coral == "P. sinensis"), cor.test(Percent_settled, nRPKM))$p.value
r_psin <- with(agmo_log.df %>% filter(., Coral == "P. sinensis"), cor.test(Percent_settled, nRPKM))$estimate

#dfav
p_dfav <- with(agmo_log.df %>% filter(., Coral == "D. favus"), cor.test(Percent_settled, nRPKM))$p.value
r_dfav <- with(agmo_log.df %>% filter(., Coral == "D. favus"), cor.test(Percent_settled, nRPKM))$estimate

#plob
p_plob <- with(agmo_log.df %>% filter(., Coral == "P. lobata"), cor.test(Percent_settled, nRPKM))$p.value
r_plob <- with(agmo_log.df %>% filter(., Coral == "P. lobata"), cor.test(Percent_settled, nRPKM))$estimate

#easp
p_easp <- with(agmo_log.df %>% filter(., Coral == "E. aspera"), cor.test(Percent_settled, nRPKM))$p.value
r_easp <- with(agmo_log.df %>% filter(., Coral == "E. aspera"), cor.test(Percent_settled, nRPKM))$estimate

# add annotations (adding stat cor manually due to formatting issues)
panel_labels <- data.frame(
  Coral = c("P. sinensis", "D. favus", "E. aspera", "P. lobata"),
  label = c("Coef = -0.172; FDR < 0.001\nR = -0.800; p < 0.001",
            "Coef = -0.002; FDR = 0.983\nR = -0.141; p = 0.421", 
            "Coef = -0.003; FDR = 0.972\nR = -0.547; p < 0.001",
            "Coef = -0.115; FDR = 0.045\nR = -0.108; p = 0.503"),
  Percent_settled = 2,   
  nRPKM = 0.28
)

panel_labels$Coral <- factor(panel_labels$Coral,
                             levels = levels(gaba.df$Coral))

p2 +
  geom_text(data = panel_labels,
            aes(x = Percent_settled, y = nRPKM, hjust = 0, label = label),
            inherit.aes = FALSE, size = 3.5)

# save
ggsave("agmo_regression_plot.svg", 
       device = "svg", 
       width = 20, 
       height = 14, 
       units = "cm", 
       dpi = 300)


