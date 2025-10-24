# Get taxonomic profile from singleM (raw reads)

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)
library(reshape2)
library(grid)
library(scales)
library(RColorBrewer)

# load functions
source("~/Documents/R/scripts/singleM_functions.R")

## import files ----

# singlem condensed otus
singlem.df <- read.table(
  "Data/Metagenomics_full/Raw/singlem/metagenome_condensed.tsv",
  sep = "\t", header = T, strip.white = T
  )

# metadata
meta.df <- read.csv(
  "Data/Metagenomics_full/Formatted/sample_metadata_filtered.csv",
  header = T, strip.white = T
  )


## format singlem data ----

# split each taxonomic rank into its own dataframe
df.l <- load_and_process_singleM_condensed_data(singlem.df)

## phylum ##

# get df
phylum.df <- m2df(df.l$taxonomy_phylum$abundances, "Phylum")

# clean colnames
colnames(phylum.df) <- gsub("_trimmed_R1", "", colnames(phylum.df))
colnames(phylum.df) <- gsub("_trimmedc000_R1", "", colnames(phylum.df))

# transpose to add metadata
rownames(phylum.df) <- phylum.df$Phylum
phylum.df <- as.matrix(select(phylum.df, -Phylum))
phylum.df <- t(phylum.df) %>% as.data.frame()

# add metadata
phylum.df$Sample_ID <- rownames(phylum.df)
phylum.df <- left_join(meta.df, phylum.df, by = "Sample_ID")

# remove phyla no longer in data (samples removed in metadata - script 1)
phylum.df <- phylum.df %>%
  dplyr::select(where(~ !is.numeric(.) || sum(.) != 0))

# filter blanks
phylum_f.df <- phylum.df %>%
  filter(phylum.df$Treatment != "negative")

# transform to long format
phylum_f.df <- phylum_f.df %>% 
  pivot_longer(-c(1:8), 
               names_to = "Phylum", 
               values_to = "Relative abundance") %>%
  as.data.frame()

## family ##

# get df
family.df <- m2df(df.l$taxonomy_family$abundances, "Family")

# clean colnames
colnames(family.df) <- gsub("_trimmed_R1", "", colnames(family.df))
colnames(family.df) <- gsub("_trimmedc000_R1", "", colnames(family.df))

# transpose to add metadata
rownames(family.df) <- family.df$Family
family.df <- as.matrix(select(family.df, -Family))
family.df <- t(family.df) %>% as.data.frame()

# add metadata
family.df$Sample_ID <- rownames(family.df)
family.df <- left_join(meta.df, family.df, by = "Sample_ID")

# remove families no longer in data (samples removed in metadata - script 1)
family.df <- family.df %>%
  dplyr::select(where(~ !is.numeric(.) || sum(.) != 0))

# filter blanks
family_f.df <- family.df %>% 
  filter(family.df$Treatment != "negative")

# make long format
family_f.df <- family_f.df %>% 
  pivot_longer(-c(1:8), 
               names_to = "Family", 
               values_to = "Relative abundance") %>% 
  as.data.frame()


## get summary statistics ----

# Function get mean and SE of for each group for a particular variable
# take a dataframe, grouping variables and response variable as input
# can take 1 or 2 groups
sum_stats <- function(df, group1, group2 = NULL, variable) {
  # format  columns (convery string to symbol)
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

## by phylum
summary_phylum <- sum_stats(phylum_f.df, 
                            "Phylum", 
                            variable = "Relative abundance")
# sort
summary_phylum <- summary_phylum %>%
  arrange(desc(mean))

# round
summary_phylum <- summary_phylum %>%
  mutate(across(where(is.numeric), ~ round(., 3)))


## by treatment and phylum
summary_treatment_p <- sum_stats(phylum_f.df, 
                                 "Treatment", 
                                 "Phylum", 
                                 variable = "Relative abundance")
# sort
summary_treatment_p <- summary_treatment_p %>%
  arrange(Treatment, desc(mean))

# round
summary_treatment_p <- summary_treatment_p %>%
  mutate(across(where(is.numeric), ~ round(., 3)))


## by family
summary_family <- sum_stats(family_f.df, 
                            "Family", 
                            variable = "Relative abundance")
# sort
summary_family <- summary_family %>%
  arrange(desc(mean))

# round
summary_family <- summary_family %>%
  mutate(across(where(is.numeric), ~ round(., 3)))


## by treatment and family
summary_treatment_f <- sum_stats(family_l.df, 
                                 "Treatment", 
                                 "Family", 
                                 variable = "Relative abundance")

# sort
summary_treatment_f <- summary_treatment_f %>%
  arrange(Treatment, desc(mean))

# round
summary_treatment_f <- summary_treatment_f %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

## save ##
write.csv(x = summary_phylum, 
          file = "phylum_relative_abundance_summary.csv", 
          row.names = F)

write.csv(x = summary_treatment_p, 
          file = "phylum_treatment_relative_abundance_summary.csv", 
          row.names = F)

write.csv(x = summary_family, 
          file = "family_relative_abundance_summary.csv", 
          row.names = F)

write.csv(x = summary_treatment_f, 
          file = "family_treatment_relative_abundance_summary.csv", 
          row.names = F)

# remove unassigned to get number of classified families
summary_family <- summary_family %>%
  filter(!str_detect(Family, "Unassigned"))

## plot top 20 phyla/family ----

## format data ##

# reload files
phylum.df <- m2df(df.l$taxonomy_phylum$abundances, "Phylum")
family.df <- m2df(df.l$taxonomy_family$abundances, "Family")

# add total column (to plot top n abundant taxa)
phylum.df$total <- rowSums(phylum.df[,-1])
family.df$total <- rowSums(family.df[,-1])

# re-order rows by total
phylum.df <- phylum.df[order(-phylum.df[,which(colnames(phylum.df) == 'total')]),]
family.df <- family.df[order(-family.df[,which(colnames(family.df) == 'total')]),]

# group families outside top 20 into 'other' category
Other <- colSums(phylum.df[-c(1:19), -1])
phylum.df <- rbind(c("Other", Other), phylum.df)

Other <- colSums(family.df[-c(1:19), -1])
family.df <- rbind(c("Other", Other), family.df)

# remove total column 
phylum.df <- select(phylum.df, !total)
family.df <- select(family.df, !total)

# subset to top n (including other)
phylum.df <- phylum.df[1:20,]
family.df <- family.df[1:20,]

# change to numeric (rbind makes it character)
tmp <- as.data.frame(sapply(phylum.df[,-1], as.numeric))
Phylum <- phylum.df$Phylum
phylum.df <- cbind(Phylum, tmp)

tmp <- as.data.frame(sapply(family.df[,-1], as.numeric))
Family <- family.df$Family
family.df <- cbind(Family, tmp)

# check 
colSums(phylum.df[,-1])
colSums(family.df[,-1])

## plot phylum ##

# clean colnames
colnames(phylum.df) <- gsub("_trimmed_R1", "", colnames(phylum.df))
colnames(phylum.df) <- gsub("_trimmedc000_R1", "", colnames(phylum.df))

# transpose to add metadata
rownames(phylum.df) <- phylum.df$Phylum
phylum.df <- as.matrix(select(phylum.df, -Phylum))
phylum.df <- t(phylum.df) %>% as.data.frame()

# add metadata
phylum.df$Sample_ID <- rownames(phylum.df)
phylum.df <- left_join(meta.df, phylum.df, by = "Sample_ID")

# remove phyla no longer in data (samples removed in metadata - script 1)
phylum.df <- phylum.df %>%
  dplyr::select(where(~ !is.numeric(.) || sum(.) != 0))

# edit phylum strings
colnames(phylum.df) <- gsub("d__.*p__", "", colnames(phylum.df))
colnames(phylum.df) <- gsub("d__", "", colnames(phylum.df))

# filter blanks
phylum_f.df <- phylum.df %>%
  filter(phylum.df$Treatment != "negative")

# transform to long format
phylum_f.df <- phylum_f.df %>% 
  pivot_longer(-c(1:8), 
               names_to = "Phylum", 
               values_to = "Relative abundance") %>% 
  as.data.frame()

# reorder treatment
phylum_f.df$Treatment <- factor(phylum_f.df$Treatment, 
                                levels = c("control", "D_2M","L_1M", "L_2M"))

# reorder phylum
phylum_f.df <- phylum_f.df %>% 
  mutate(Phylum = fct_relevel(Phylum,
                              "Bacteria;Unassigned", 
                              "Archaea;Unassigned", 
                              "Other",
                              after = Inf))

# get facet labels
facet_names <- c(
  # `negative` = "Blank",
  `control` = "Ctrl",
  `L_1M` = "Light 1M",
  `L_2M` = "Light 2M",
  `D_2M` = "Dark 2M"
)

# get colour palette
cols <- colorRampPalette(brewer.pal(12, "Paired"))(20)
show_col(cols)

# plot stacked bar
p <- ggplot(phylum_f.df, aes(fill=Phylum, colour=Phylum, 
                             y=`Relative abundance`, 
                             x=Sample_ID)) + 
  geom_bar(position="fill", stat="identity", width=1) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  xlab("Biofilm Sample") +
  ylab("Relative Abundance") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.position='right') 
p + facet_grid(~Treatment, 
               scales = 'free_x', 
               space = 'free', 
               labeller = as_labeller(facet_names)) + 
  theme(panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1))

# save
ggsave("singlem_rel_abun_phylum_bar_t20.svg", 
       device = "svg", 
       width = 26, 
       height = 15, 
       dpi = 300, 
       units = "cm")

## plot family ##

# clean colnames
colnames(family.df) <- gsub("_trimmed_R1", "", colnames(family.df))
colnames(family.df) <- gsub("_trimmedc000_R1", "", colnames(family.df))

# transpose to add metadata
rownames(family.df) <- family.df$Family
family.df <- as.matrix(select(family.df, -Family))
family.df <- t(family.df) %>% as.data.frame()

# add metadata
family.df$Sample_ID <- rownames(family.df)
family.df <- left_join(meta.df, family.df, by = "Sample_ID")

# remove families no longer in data (samples removed in metadata - script 1)
family.df <- family.df %>%
  dplyr::select(where(~ !is.numeric(.) || sum(.) != 0))

# edit family strings
colnames(family.df) <- gsub("d__.*p__.*c__.*o__", "", colnames(family.df))
colnames(family.df) <- gsub("d__Bacteria;Unassigned.*", "Bacteria;Unassigned", colnames(family.df))
colnames(family.df) <- gsub(".*Alphaproteobacteria;Unassigned.*", "Alphaproteobacteria;Unassigned", colnames(family.df))
colnames(family.df) <- gsub(".*Gammaproteobacteria;Unassigned.*", "Gammaproteobacteria;Unassigned", colnames(family.df))
colnames(family.df) <- gsub(".*;f__", "", colnames(family.df))
colnames(family.df) <- gsub("MH13", "Rhizobiales;MH13", colnames(family.df))
colnames(family.df) <- gsub("HXMU1428-3", "Rhizobiales;HXMU1428-3", colnames(family.df))

# filter blanks
family_f.df <- family.df %>% 
  filter(family.df$Treatment != "negative")

# make long format
family_f.df <- family_f.df %>% 
  pivot_longer(-c(1:8), 
               names_to = "Family", 
               values_to = "Relative abundance") %>% 
  as.data.frame()

# reorder treatment
family_f.df$Treatment <- factor(family_f.df$Treatment, 
                              levels = c("control","D_2M","L_1M", "L_2M"))
# reorder family 
family_f.df <- family_f.df %>% 
  mutate(Family = fct_relevel(Family, 
                              "Bacteria;Unassigned",
                              "Other", 
                              after = Inf))

# plot stacked bar
p <- ggplot(family_f.df, aes(fill=Family, colour=Family, 
                             y=`Relative abundance`, 
                             x=Sample_ID)) + 
  geom_bar(position="fill", stat="identity", width=1) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  xlab("Biofilm Sample") +
  ylab("Relative Abundance") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.position='right') 
p + facet_grid(~Treatment, 
               scales = 'free_x', 
               space = 'free', 
               labeller = as_labeller(facet_names)) + 
  theme(panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1))

ggsave("singlem_rel_abun_family_t20_bar.svg", 
       device = "svg", 
       width = 28.4, 
       height = 15, 
       dpi = 300, 
       units = "cm")

