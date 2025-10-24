# Get gene centric alpha diversity - running rstudio on bunya server

setwd("/scratch/project/micro_inducers")

library(tidyverse)
library(vegan)
library(viridis)
library(scales)
library(car)
library(multcomp)

## load and format data ----

genes.df <-read.table(
  "analysis/ECT03-1_biofilm_MG_2022/08_annotate/cdhit/rpkm/combined_rpkm_results/all_samples_cdhit100_rpkm_per_gene_normalised.tsv", 
  header = T, sep = "\t", strip.white = T, fill = T
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

# get rownames and transpose
rownames(genes.df) <- genes.df$Gene_ID
genes.df <- t(genes.df[,-1])

# get gene richness
rich.df <- apply(genes.df[,-1] >0, 1, sum) %>% 
  as.data.frame() %>%
  dplyr::rename(richness = 1) 
rich.df$Sample_ID <- rownames(rich.df)

# get shannons diversity
shan.df <- diversity(genes.df, index = "shannon") %>% 
  as.data.frame %>%
  dplyr::rename(shannon = 1) 
shan.df$Sample_ID <- rownames(shan.df)

# join with meta
rich.df <-left_join(meta.df, rich.df, by = "Sample_ID")
shan.df <- left_join(meta.df, shan.df, by = "Sample_ID")

# reorder treatment
rich.df$Treatment <- factor(
  rich.df$Treatment, levels = c("control", "D_2M", "L_1M", "L_2M")
)

shan.df$Treatment <- factor(
  shan.df$Treatment, levels = c("control", "D_2M", "L_1M", "L_2M")
)

## calculate significance values ----

# transform data
rich.df$rich_fourth <- (rich.df$richness ^ (1/4))
shan.df$shan_fourth <- (shan.df$shannon ^ (1/4))

# check assumptions 
#richness
residualPlots(lm(rich.df$richness~rich.df$Treatment))
residualPlots(lm(rich.df$rich_fourth~rich.df$Treatment)) 

#shannon
residualPlots(lm(shan.df$shannon~shan.df$Treatment))
residualPlots(lm(shan.df$shan_fourth~shan.df$Treatment)) 

# get lm
#richness
rich.lm <- lm(rich_fourth ~ Treatment, rich.df)
plot(rich.lm)

# summarise linear model. Note: treatment names are retained when factor levels are not re-ordered
summary(rich.lm) 
anova(rich.lm) # for anova

#shannon
shan.lm <- lm(shan_fourth ~ Treatment, shan.df)
plot(shan.lm)

# summarise
summary(shan.lm)
anova(shan.lm) 

#post hoc test
# richness
rich.glht <- summary(glht(rich.lm, linfct = mcp(Treatment = "Tukey")), test = adjusted("bonferroni"))
rich.groups <- cld(rich.glht)
rich.groups <- fortify(rich.groups)
colnames(rich.groups) <- c("Treatment", "letters")

ymax <- tapply(rich.df$richness, rich.df$Treatment, max)
rich.groups$Ymax <- ymax # add to plot below

# shannon
shan.glht <- summary(glht(shan.lm, linfct = mcp(Treatment = "Tukey")), test = adjusted("bonferroni"))
shan.groups <- cld(shan.glht)
shan.groups <- fortify(shan.groups)
colnames(shan.groups) <- c("Treatment", "letters")

ymax <- tapply(shan.df$shannon, shan.df$Treatment, max)
shan.groups$Ymax <- ymax # add to plot below

## save dataframes
saveRDS(rich.df, "gene_richness_df.rds")
saveRDS(rich.groups, "genes_richness_groups.rds")
saveRDS(shan.df, "gene_shannon_df.rds")
saveRDS(shan.groups, "gene_shannon_groups.rds")

## plot diversity ----

# get colours 
clrs <- viridis(4, direction = 1, option = "turbo")
clrs <- clrs[c(1:3)]
clrs <- c("#808080", clrs)
show_col(clrs)

# plot richness
p <- ggplot(rich.df, aes(x=Treatment, y=richness, fill=Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  ylab("Gene richness") +
  scale_fill_manual(values = clrs) +
  scale_x_discrete(labels = c("Control", "Dark 2M", "Light 1M", "Light 2M")) +
  guides(fill = guide_legend(title = "Treatment")) +
  theme_bw() +
  theme(axis.text.x=element_text(size=11, angle = 45, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 13),
        legend.position = "none") 
p

# to add significance values
p + geom_text(data=rich.groups, 
              aes(x = Treatment, y = Ymax + 300000, label = letters),
              vjust=0, size = 4)
# save
ggsave("gene_richness_boxplot.svg", device = "svg", width = 15, height = 12, 
       units = "cm", dpi = 300)

# plot shannon
p <- ggplot(shan.df, aes(x=Treatment, y=shannon, fill=Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  ylim(5, NA) +
  ylab("Shannon Diversity Index") +
  scale_fill_manual(values = clrs) +
  scale_x_discrete(labels = c("Control", "Dark 2M", "Light 1M", "Light 2M")) +
  guides(fill = guide_legend(title = "Treatment")) +
  theme_bw() +
  theme(axis.text.x=element_text(size=11, angle = 45, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 13),
        legend.position = "none") 
p

# to add significance values
p + geom_text(data=shan.groups, 
              aes(x = Treatment, y = Ymax+0.75, label = letters),
              vjust=0, size = 4)
# save
ggsave("gene_shannon_boxplot.svg", device = "svg", width = 15, height = 12, 
       units = "cm", dpi = 300)



