## Get beta diversity of taxonomy using MAGs and singleM results

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)
library(vegan)
library(viridis)
library(scales)

# load singlem functions
source("~/Documents/R/example_scripts/singleM_functions.R")

## import and format data ----

## singleM
singlem.df <- read.table(
  "Data/Metagenomics_full/Raw/singlem/metagenome_condensed.tsv",
  sep = "\t", header = T, strip.white = T
)

# split each taxonomic rank into its own dataframe
df.l <- load_and_process_singleM_condensed_data(singlem.df)

# convert to species table
species.df <- m2df(df.l$taxonomy_species$abundances, "Species")

# clean colnames
colnames(species.df) <- gsub("_trimmed_R1", "", colnames(species.df))
colnames(species.df) <- gsub("_trimmedc000_R1", "", colnames(species.df))

### mags
mag.df <- read.table(
  "Data/Metagenomics_full/Formatted/dereplicated_mags_summary_data_filtered.tsv",
  sep = '\t', header = T, strip.white = T
)

# remove unmapped
mag.df <- mag.df %>%
  filter(MAG_ID != "unmapped")

# get mag metadata
meta_mags.df <- mag.df %>%
  dplyr::select(!starts_with("SE")) 

# subset to mag_id and samples columns
mag.df <- mag.df %>%
  dplyr::select(MAG_ID, starts_with("SE")) 

## sample metadata
meta.df <- read.csv(
  "Data/Metagenomics_full/Formatted/sample_metadata_filtered.csv",
  header = T, strip.white = T
)

# remove blanks 
meta.df <- meta.df %>%
  filter(Treatment != "negative")

mag.df <- mag.df %>%
  dplyr::select(MAG_ID, all_of(meta.df$Sample_ID))

species.df <- species.df %>%
  dplyr::select(Species, all_of(meta.df$Sample_ID))

## transpose to make matrix
#mags
row.names(mag.df) <- mag.df$MAG_ID
mag.mat <- mag.df %>%
  dplyr::select(!MAG_ID) %>%
  t()

#singlem
row.names(species.df) <- species.df$Species
species.mat <- species.df %>%
  dplyr::select(!Species) %>%
  t()


## Run NMDS ----

# function for nmds plot 
plot_nmds <- function(df, distance, trymax, autotransform, k, plot, metadata) {
  # run nmds
  dist.nmds <- metaMDS(df, 
                       distance = distance, 
                       trymax = trymax, 
                       autotransform = autotransform, 
                       k = k, 
                       plot = plot)
  # Extract scores from metaMDS and return as data matrix
  nmds.scores <- scores(dist.nmds)
  nmds.scores.site <- as.data.frame(nmds.scores$sites)
  # add metadata
  nmds.scores.site$Sample_ID <- row.names(nmds.scores.site)
  nmds.df <- left_join(nmds.scores.site, metadata, by = "Sample_ID")
  return(nmds.df)
}

# run function over matrix
# mag
mag.nmds <- plot_nmds(df = mag.mat, 
                      distance = "bray",
                      trymax = 500,
                      autotransform = T,
                      k = 2,
                      plot = T,
                      metadata = meta.df)

# species
spec.nmds <- plot_nmds(df = species.mat, 
                      distance = "bray",
                      trymax = 500,
                      autotransform = T,
                      k = 2,
                      plot = T,
                      metadata = meta.df)

## Format NMDS and plot ----

# mags
# reorder treatment
mag.nmds$Treatment <- factor(
  mag.nmds$Treatment, levels = c("control", "D_2M", "L_1M", "L_2M")
)

# treatment colours
clrs <- viridis(4, direction = 1, option = "turbo")
clrs <- clrs[c(1:3)]
clrs <- c("#808080", clrs)
show_col(clrs)

# get labels
x.labs <- c("Control", "Dark 2M", "Light 1M", "Light 2M")

# plot
p <- ggplot(mag.nmds, aes(x=NMDS1, y=NMDS2, 
                          colour = Treatment,
                          shape = Settlement_2,
                          label = Sample_ID)) +
  geom_point(size = 4, alpha = 0.7) + 
  #geom_text(hjust=1.3, vjust=0) +
  scale_colour_manual(values = clrs,
                      labels = c("Control", "Dark 2M", "Light 1M", "Light 2M")) +
  labs(colour = "Treatment",
       shape = "Settlement") +
  guides(shape = guide_legend(order = 2), 
         colour = guide_legend(order = 1)) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14)) 
p

# save
ggsave("nmds_mags.svg",
       device = "svg", 
       width = 20, 
       height = 15, 
       dpi = 300, 
       units = "cm")

# species
# reorder treatment
spec.nmds$Treatment <- factor(
  spec.nmds$Treatment, levels = c("control", "D_2M", "L_1M", "L_2M")
)

# plot
p <- ggplot(spec.nmds, aes(x=NMDS1, y=NMDS2, 
                           colour = Treatment,
                           shape = Settlement_2,
                           label = Sample_ID)) +
  geom_point(size = 4, alpha = 0.7) + 
  #geom_text(hjust=1.3, vjust=0) +
  scale_colour_manual(values = clrs,
                      labels = c("Control", "Dark 2M", "Light 1M", "Light 2M")) +
  labs(colour = "Treatment",
       shape = "Settlement") +
  guides(shape = guide_legend(order = 2), 
         colour = guide_legend(order = 1)) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14)) 
p

# save
ggsave("nmds_species.svg",
       device = "svg", 
       width = 20, 
       height = 15, 
       dpi = 300, 
       units = "cm")

## Run PERMANOVA ----

## create distance matrix 
# standardise data
mag_dw.mat <- wisconsin(mag.mat)
spec_dw.mat <- wisconsin(species.mat)

# get distance
mag.dist <- vegdist(mag_dw.mat)
spec.dist <- vegdist(spec_dw.mat)

# run PERMANOVA
# mag
mag.perm <- adonis2(mag.dist~Treatment, 
                    meta.df, 
                    permutations = 999)
mag.perm

# species
spec.perm <- adonis2(spec.dist~Treatment, 
                     meta.df, 
                     permutations = 999)
spec.perm
  
# save results
saveRDS(mag.perm, "Results/metagenomics/permanova/mag_permanova_treatment.rds")
saveRDS(spec.perm, "Results/metagenomics/permanova/species_permanova_treatment.rds")







