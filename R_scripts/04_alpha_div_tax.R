# Get alpha diversity of taxonomy using MAGs and singleM results

setwd("~/Documents/R/Microbial_inducers")

library(tidyverse)
library(vegan)
library(viridis)
library(scales)
library(car)
library(multcomp)

# load functions
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

## mags
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

## transpose
#mags
row.names(mag.df) <- mag.df$MAG_ID
mag.df <- mag.df %>%
  dplyr::select(!MAG_ID) %>%
  t()

#singlem
row.names(species.df) <- species.df$Species
species.df <- species.df %>%
  dplyr::select(!Species) %>%
  t()


## calculate and plot species richness ----

# calculate richness
# mags
rich.m <- apply(mag.df >0, 1, sum) %>% 
  as.data.frame() %>%
  rename("mag_richness" = ".") %>%
  mutate(Sample_ID = rownames(.))

# species
rich.s <- apply(species.df >0, 1, sum) %>%
  as.data.frame() %>%
  rename("species_richness" = ".") %>% 
  mutate(Sample_ID = rownames(.))

# combine with meta
richness.df <- left_join(meta.df, rich.m)
richness.df <- left_join(richness.df, rich.s)

# remove blanks
richness.df <- richness.df %>%
  filter(Treatment != "negative")

# reorder treatment
richness.df$Treatment <- factor(
  richness.df$Treatment, levels = c("control", "D_2M", "L_1M", "L_2M")
)

# calculate significance values 
# mags
# transform data
richness.df$fourth_mag <- (richness.df$mag_richness ^ (1/4))

# check assumptions 
residualPlots(lm(richness.df$mag_richness~richness.df$Treatment)) 
residualPlots(lm(richness.df$fourth_mag~richness.df$Treatment)) # better to use transformed data as less variance between treatments

# get lm
mag_rich.lm <- lm(fourth_mag ~ Treatment, richness.df)
plot(mag_rich.lm)

summary(mag_rich.lm) # summarise linear model. Note: treatment names are retained when factor levels are not re-ordered
anova(mag_rich.lm) # for anova

#post hoc test
mag_rich.glht <- summary(glht(mag_rich.lm, linfct = mcp(Treatment = "Tukey")), test = adjusted("bonferroni"))
mag_rich.groups <- cld(mag_rich.glht)
mag_rich.groups <- fortify(mag_rich.groups)
colnames(mag_rich.groups) <- c("Treatment", "letters")

ymax <- tapply(richness.df$mag_richness, richness.df$Treatment, max)
mag_rich.groups$Ymax <- ymax # add to plot 

# singlem
# transform data
richness.df$fourth_spe <- (richness.df$species_richness ^ (1/4))

# check assumptions 
residualPlots(lm(richness.df$species_richness~richness.df$Treatment)) 
residualPlots(lm(richness.df$fourth_spe~richness.df$Treatment)) # better to use transformed data as less variance between treatments

# get lm
spec_rich.lm <- lm(fourth_spe ~ Treatment, richness.df)
plot(spec_rich.lm)

summary(spec_rich.lm) # summarise linear model. Note: treatment names are retained when factor levels are not re-ordered
anova(spec_rich.lm) # for anova

#post hoc test
spec_rich.glht <- summary(glht(spec_rich.lm, linfct = mcp(Treatment = "Tukey")), test = adjusted("bonferroni"))
spec_rich.groups <- cld(spec_rich.glht)
spec_rich.groups <- fortify(spec_rich.groups)
colnames(spec_rich.groups) <- c("Treatment", "letters")

ymax <- tapply(richness.df$species_richness, richness.df$Treatment, max)
spec_rich.groups$Ymax <- ymax # add to plot 


## plot boxplot
# treatment colours
clrs <- viridis(4, direction = 1, option = "turbo")
clrs <- clrs[c(1:3)]
clrs <- c("#808080", clrs)
show_col(clrs)

# get labels
x.labs <- c("Control", "Dark 2M", "Light 1M", "Light 2M")

# mags
p <- ggplot(richness.df, aes(x=Treatment, y=mag_richness, fill=Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  ylab("MAG richness") +
  scale_fill_manual(values = clrs,
                    labels = x.labs) +
  scale_x_discrete(labels = x.labs) +
  guides(fill = guide_legend(title = "Treatment")) +
  theme_bw() +
  theme(axis.text.x=element_text(size=11),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 13),
        legend.position = "none") 
p

# to add significance values
p + geom_text(data=mag_rich.groups, 
              aes(x = Treatment, y = Ymax+25, label = letters),
              vjust=0, size = 5)

# save
ggsave("mag_richness.svg", device = "svg", width = 15, height = 12, 
       units = "cm", dpi = 300)

# singlem
p <- ggplot(richness.df, aes(x=Treatment, y=species_richness, fill=Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  ylab("Species richness") +
  scale_fill_manual(values = clrs,
                    labels = x.labs) +
  scale_x_discrete(labels = x.labs) +
  guides(fill = guide_legend(title = "Treatment")) +
  theme_bw() +
  theme(axis.text=element_text(size=11),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 13),
        legend.position = "none") 
p
# to add significance values - see below
p + geom_text(data=spec_rich.groups, 
              aes(x = Treatment, y = Ymax+75, label = letters),
              vjust=0, size = 5)
# save
ggsave("species_richness.svg", device = "svg", width = 15, height = 12, 
       units = "cm", dpi = 300)

## get group averages 

# function for std error
std <- function(x) sd(x)/sqrt(length(x))

# mags
mag_rich_mean <- richness.df %>%
  group_by(Treatment) %>%
  dplyr::summarise(Mean = mean(mag_richness, na.rm=TRUE),
                   SD = sd(mag_richness, na.rm=TRUE),
                   SE = std(mag_richness))

# species
spec_rich_mean <- richness.df %>%
  group_by(Treatment) %>%
  dplyr::summarise(Mean = mean(species_richness, na.rm=TRUE),
                   SD = sd(species_richness, na.rm=TRUE),
                   SE = std(species_richness))

## calculate and plot shannons diversity ----

# mags
shan.m <- diversity(mag.df, index = "shannon") %>%
  as.data.frame() %>%
  rename("mag_shannon" = ".") %>%
  mutate(Sample_ID = rownames(.))

# singlem
shan.s <- diversity(species.df, index = "shannon") %>%
  as.data.frame() %>%
  rename("spec_shannon" = ".") %>%
  mutate(Sample_ID = rownames(.))

# combine with meta
shannon.df <- left_join(meta.df, shan.m)
shannon.df <- left_join(shannon.df, shan.s)

# remove blanks
shannon.df <- shannon.df %>%
  filter(Treatment != "negative")

# reorder treatment
shannon.df$Treatment <- factor(
  shannon.df$Treatment, levels = c("control", "D_2M", "L_1M", "L_2M")
)

# calculate significance values 
# mags
# transform data
shannon.df$fourth_mag <- (shannon.df$mag_shannon ^ (1/4))

# check assumptions 
residualPlots(lm(shannon.df$mag_shannon~shannon.df$Treatment)) 
residualPlots(lm(shannon.df$fourth_mag~shannon.df$Treatment))

# get lm
mag_shan.lm <- lm(fourth_mag ~ Treatment, shannon.df)
plot(mag_shan.lm)

# summarise model
summary(mag_shan.lm)
anova(mag_shan.lm)

#post hoc test
mag_shan.glht <- summary(glht(mag_shan.lm, linfct = mcp(Treatment = "Tukey")), test = adjusted("bonferroni"))
mag_shan.groups <- cld(mag_shan.glht)
mag_shan.groups <- fortify(mag_shan.groups)
colnames(mag_shan.groups) <- c("Treatment", "letters")

ymax <- tapply(shannon.df$mag_shannon, shannon.df$Treatment, max)
mag_shan.groups$Ymax <- ymax # add to plot 

# species
# transform data
shannon.df$fourth_spec <- (shannon.df$spec_shannon ^ (1/4))

# check assumptions 
residualPlots(lm(shannon.df$spec_shannon~shannon.df$Treatment)) # untransformed more suitable
residualPlots(lm(shannon.df$fourth_spec~shannon.df$Treatment))

# get lm
spec_shan.lm <- lm(spec_shannon ~ Treatment, shannon.df)
plot(spec_shan.lm)

# summarise model
summary(spec_shan.lm)
anova(spec_shan.lm)

#post hoc test
spec_shan.glht <- summary(glht(spec_shan.lm, linfct = mcp(Treatment = "Tukey")), test = adjusted("bonferroni"))
spec_shan.groups <- cld(spec_shan.glht)
spec_shan.groups <- fortify(spec_shan.groups)
colnames(spec_shan.groups) <- c("Treatment", "letters")

ymax <- tapply(shannon.df$spec_shannon, shannon.df$Treatment, max)
spec_shan.groups$Ymax <- ymax # add to plot 

## plot boxplot
# treatment colours
clrs <- viridis(4, direction = 1, option = "turbo")
clrs <- clrs[c(1:3)]
clrs <- c("#808080", clrs)
show_col(clrs)

# get labels
x.labs <- c("Control", "Dark 2M", "Light 1M", "Light 2M")

# mags
p <- ggplot(shannon.df, aes(x=Treatment, y=mag_shannon, fill=Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  ylab("Shannon Diversity Index (MAGs)") +
  scale_fill_manual(values = clrs,
                    labels = x.labs) +
  scale_x_discrete(labels = x.labs) +
  guides(fill = guide_legend(title = "Treatment")) +
  theme_bw() +
  theme(axis.text.x=element_text(size=11),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 13),
        legend.position = "none") 
p

# to add significance values - see below
p + geom_text(data=mag_shan.groups, 
              aes(x = Treatment, y = Ymax+0.2, label = letters),
              vjust=0, size = 5)

# save
ggsave("mag_shannon.svg", device = "svg", width = 15, height = 12, 
       units = "cm", dpi = 300)

# singlem
p <- ggplot(shannon.df, aes(x=Treatment, y=spec_shannon, fill=Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  ylab("Shannon Diversity Index (SingleM)") +
  scale_fill_manual(values = clrs,
                    labels = x.labs) +
  scale_x_discrete(labels = x.labs) +
  guides(fill = guide_legend(title = "Treatment")) +
  theme_bw() +
  theme(axis.text=element_text(size=11),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 13),
        legend.position = "none") 
p
# to add significance values - see below
p + geom_text(data=spec_shan.groups, 
              aes(x = Treatment, y = Ymax+0.2, label = letters),
              vjust=0, size = 5)
# save
ggsave("species_shannon.svg", device = "svg", width = 15, height = 12, 
       units = "cm", dpi = 300)


## get group averages 
# mags
mag_shan_mean <- shannon.df %>%
  group_by(Treatment) %>%
  dplyr::summarise(Mean = mean(mag_shannon, na.rm=TRUE),
                   SD = sd(mag_shannon, na.rm=TRUE),
                   SE = std(mag_shannon))

# species
spec_shan_mean <- shannon.df %>%
  group_by(Treatment) %>%
  dplyr::summarise(Mean = mean(spec_shannon, na.rm=TRUE),
                   SD = sd(spec_shannon, na.rm=TRUE),
                   SE = std(spec_shannon))





