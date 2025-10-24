# Get gene centric beta diversity - running rstudio on bunya server

setwd("/scratch/project/micro_inducers")

library(tidyverse)
library(vegan)
library(viridis)
library(scales)
#library(parallelDist)
library(car)
library(multcomp)


## load and format data ----

# load genes df
genes.df <- read.table(
  "/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/08_annotate/cdhit/rpkm/combined_rpkm_results/all_samples_cdhit100_rpkm_per_gene_normalised.tsv", 
  header = T,
  sep = "\t",
  strip.white = T
)

# load metadata
meta.df <- read.csv(
  "data/metadata_files/sample_metadata_filtered.csv",
  header = T, strip.white = T
)

# remove blanks
meta.df <- meta.df %>%
  filter(Treatment != "negative")

# subset samples to meta
genes.df <- genes.df %>%
  dplyr::select(Gene_ID, all_of(meta.df$Sample_ID))

# remove any genes no longer in df
genes.df <- genes.df[rowSums(genes.df[,-1]) > 0,]

## run NMDS ----

# create data matrix
mat <- genes.df # change to necessary dataframe
rownames(mat) <- genes.df$Gene_ID
mat <- t(mat[,-1])

# transform data (high variability and skewness)
mat.log <- log1p(mat)
#mat.sqrt <- sqrt(mat)

# run nmds
bray.nmds <- metaMDS(mat.log, 
                     distance = "bray", 
                     trymax = 1000, 
                     autotransform = F, 
                     k = 2, 
                     plot = T)

# save nmds
saveRDS(bray.nmds, "/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/11_beta_diversity/bray_nmds_log.rds")
bray.nmds <- readRDS("/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/11_beta_diversity/bray_nmds_log.rds")

## format NMDS and plot ----

# extract scores from metaMDS and return as data matrix
nmds.scores <- scores(bray.nmds)
nmds.scores.site <- as.data.frame(nmds.scores$sites)

# add metadata
nmds.scores.site$Sample_ID <- row.names(nmds.scores.site)
nmds.df <- left_join(meta.df, nmds.scores.site, by = "Sample_ID")

# remove outlier (ctrl with no reads)
nmds.df <- nmds.df %>%
  filter(Sample_ID != "SE2469")

## plot nmds
# get colours 
clrs <- viridis(4, direction = 1, option = "turbo")
clrs <- clrs[c(1:3)]
clrs <- c("#808080", clrs)

# relevel settlement
nmds.df <- nmds.df %>%
  mutate(Settlement = fct_relevel(Settlement, "None", "Low", "High"))

# colour = treatment, shape = settlement
p <- ggplot(nmds.df, aes(x=NMDS1, y=NMDS2, 
                         colour = Treatment,
                         shape = Settlement,
                         label = Sample_ID)) +
  geom_point(size = 4, alpha = 0.7) + 
  # geom_text(hjust=1.3, vjust=0) +
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
ggsave(
  "/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/11_beta_diversity/nmds_genes_log_treatment_settlement.svg",
  device = "svg", width = 20, height = 15, 
  dpi = 300, units = "cm")

# colour = treatment, shape = tank
p <- ggplot(nmds.df, aes(x=NMDS1, y=NMDS2, 
                         colour = Treatment,
                         fill = Treatment,
                         shape = Tank)) +
  geom_point(size = 4, alpha = 0.7) + 
  scale_fill_manual(values = clrs,
                    labels = c("Control", "Dark 2M", "Light 1M", "Light 2M")) +
  scale_colour_manual(values = clrs,
                      labels = c("Control", "Dark 2M", "Light 1M", "Light 2M")) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("Control", "Tank 1", "Tank 2", "Tank 3")) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14))
p

# save
ggsave(
  "/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/11_beta_diversity/nmds_genes_log_treatment_tank.svg",
  device = "svg", width = 20, height = 15, 
  dpi = 300, units = "cm")

# colour = settlement, shape = treatment
p <- ggplot(nmds.df, aes(x=NMDS1, y=NMDS2, 
                         colour = Percent_settled,
                         fill = Percent_settled,
                         shape = Treatment)) +
  geom_point(size = 4, alpha = 0.7) + 
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("Control", "Dark 2M", "Light 1M", "Light 2M")) +
  scale_colour_viridis(discrete = F, option = "cividis", direction = 1) +
  scale_fill_viridis(discrete = F, option = "cividis", direction = 1) +
  labs(colour = "Settlement (%)",
       fill = "Settlement (%)",
       shape = "Treatment") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14))
p

# save
ggsave(
  "/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/11_beta_diversity/nmds_genes_sqrt_settlement_treatment.svg",
  device = "svg", width = 20, height = 13, 
  dpi = 300, units = "cm")


## run PERMANOVA and PERMDISP ----

# remove controls to only assess differences among treatments
meta_f.df <- meta.df %>%
  filter(Treatment != "control")
genes.df <- genes.df %>%
  dplyr::select(Gene_ID, all_of(meta_f.df$Sample_ID))

# create distance matrix
# get data matrix
mat <- genes.df
rownames(mat) <- genes.df$Gene_ID
mat <- t(mat[,-1])

# transform
mat.log <- log1p(mat)

# get distance matrix
bray.dist.log <- vegdist(mat.log)

# save
saveRDS(bray.dist.log, "/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/11_beta_diversity/bray_dist_cd100_log.rds")
bray.dist.log <- readRDS("/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/11_beta_diversity/bray_dist_cd100_log.rds")

# join df with meta 
perm_log.df <- as.data.frame(mat.log)
perm_log.df$Sample_ID <- row.names(mat.log)

#perm.df <- left_join(meta_f, perm.df, by = "Sample_ID") #join not working due to large file size, using cbind instead

# ensure order matches
meta_aligned <- meta_f.df[match(perm_log.df$Sample_ID, meta_f.df$Sample_ID), ]
identical(meta_aligned$Sample_ID, perm_log.df$Sample_ID)

# cbind
perm_log.df <- cbind(meta_aligned, dplyr::select(perm_log.df, !Sample_ID))

# run permanova
perm_results <- adonis2(
  bray.dist.log~Treatment*Settlement*Tank*Coral,
  perm_log.df,
  permutations = 9999
)
perm_results

# save
saveRDS(perm_results, "/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/11_beta_diversity/permanova_results_log.rds")

# run permdisp
disp_trmt.log <- betadisper(bray.dist.log, group = meta_f.df$Treatment)
anova(disp_trmt.log)
pt_trmt.log <- permutest(disp_trmt.log)

disp_set.log <- betadisper(bray.dist.log, group = meta_f.df$Settlement)
anova(disp_set.log)
pt_set.log <- permutest(disp_set.log)

disp_tank.log <- betadisper(bray.dist.log, group = meta_f.df$Tank)
anova(disp_tank.log)
pt_tank.log <- permutest(disp_tank.log)

disp_cor.log <- betadisper(bray.dist.log, group = meta_f.df$Coral)
anova(disp_cor.log)
pt_cor.log <- permutest(disp_cor.log)

# save
permdisp_results_log.l <- list(
  "Treatment" = pt_trmt.log,
  "Settlement" = pt_set.log,
  "Tank" = pt_tank.log,
  "Coral" = pt_cor.log
)
saveRDS(permdisp_results_log.l, "/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/11_beta_diversity/permdisp_results_log.rds")


