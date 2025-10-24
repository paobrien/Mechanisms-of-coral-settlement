## Get taxonomic profile of MAGs

setwd("~/Documents/R/Microbial_inducers")

library(tidyverse)
library(viridis)
library(scales)

## import MAG data ----

# mags df
mag.df <- read.table(
  "Data/Metagenomics_full/Formatted/dereplicated_mags_summary_data_filtered.tsv",
  sep = '\t', header = T, strip.white = T
  )

# remove unmapped
mag.df <- mag.df %>%
  filter(MAG_ID != "unmapped")

# metadata
meta.df <- read.csv(
  "Data/Metagenomics_full/Formatted/sample_metadata_filtered.csv",
  header = T, strip.white = T
  )


## format and plot ----

# change NAs (no mapping) to 0
mag.df[is.na(mag.df)] <- 0

# transpose so mags are columns
rownames(mag.df) <- mag.df$MAG_ID
mag_t.df <- mag.df %>%
  select(starts_with("SE")) %>%
  t() %>% as.data.frame()

# join with meta
mag_t.df$Sample_ID <- row.names(mag_t.df)
mag_t.df <- full_join(meta.df, mag_t.df, by = "Sample_ID")

# filter blanks
mag_f.df <- mag_t.df %>%
  filter(Treatment != "negative")

# make long format
mag_f.df <- mag_f.df %>% 
  pivot_longer(
    cols = starts_with("SE", ignore.case = FALSE), 
    names_to = "MAG") %>% 
  as.data.frame()

# add taxonomy
mag_f.df <- left_join(
  mag_f.df, select(mag.df, MAG_ID, Phylum:Species), join_by("MAG" == "MAG_ID"))

# reorder treatment
mag_f.df$Treatment <- factor(
  mag_f.df$Treatment, levels = c("control", "D_2M","L_1M", "L_2M")
  )

# group by phylum
p.df <- mag_f.df %>% 
  group_by(Sample_ID, Phylum, Treatment, Tank) %>%
  summarise(total = sum(value, na.rm = T))

# edit phylum strings
p.df$Phylum <- gsub("p__", "", p.df$Phylum)

# group by family
f.df <- mag_f.df %>% 
  group_by(Sample_ID, Family, Treatment, Tank) %>% 
  summarise(total = sum(value, na.rm = T))

# edit family strings
f.df$Family <- gsub("f__", "", f.df$Family)

# get range
range(p.df$total)
p.r_val <- c(0, 30, 62)

range(f.df$total)
f.r_val <- c(0, 15, 33)

# get colours
heat_col <- viridis(n = 9, 
                    alpha = 0.9, 
                    begin = 0, 
                    end = 0.95, 
                    option = "turbo", 
                    direction = 1)
show_col(heat_col)

# get facet labels
facet_names <- c(
  `control` = "Ctrl",
  `L_1M` = "Light 1M",
  `L_2M` = "Light 2M",
  `D_2M` = "Dark 2M"
)

## plot 
# phyla
p <- ggplot(p.df, aes(x = Sample_ID, y = Phylum, fill = total)) +
  geom_raster() +
  scale_fill_gradientn(name = "Relative Abundance", 
                       colours = heat_col) +#,
                       #values = rescale(p.r_val)) + 
  xlab("Biofilm Sample") +
  ylab("Phylum") +
  theme(#axis.text.x = element_text(size = 12, angle = 90),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.ticks.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12)) 
# facet grid (treatment)
p + facet_grid(~Treatment, 
               scales = 'free_x', 
               space = 'free', 
               labeller = as_labeller(facet_names)) + 
  theme(panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(color = "black", linewidth = 1))

## save
ggsave("mag_rel_abundnace_phylum.svg", device = "svg", width = 40, height = 25, 
       dpi = 300, units = "cm", limitsize = FALSE)


## family
p <- ggplot(f.df, aes(x = Sample_ID, y = Family, fill = total)) +
  geom_raster() +
  scale_fill_gradientn(name = "Relative Abundance", 
                       colours = heat_col) +#,
  #values = rescale(p.r_val)) + 
  xlab("Biofilm Sample") +
  ylab("Phylum") +
  theme(#axis.text.x = element_text(size = 12, angle = 90),
    axis.text.x = element_blank(),
    axis.title.x = element_text(size = 15),
    axis.ticks.x = element_blank(),
    #axis.text.y = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12)) 
# facet grid (treatment)
p + facet_grid(~Treatment, 
               scales = 'free_x', 
               space = 'free', 
               labeller = as_labeller(facet_names)) + 
  theme(panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(color = "black", linewidth = 1))

# save
ggsave("mag_rel_abundnace_family.svg", device = "svg", width = 40, height = 40, 
       dpi = 300, units = "cm", limitsize = FALSE)

## summarise data ----

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
summary_phylum <- sum_stats(p.df, 
                            "Phylum", 
                            variable = "total")
# sort
summary_phylum <- summary_phylum %>%
  arrange(desc(mean))

# round
summary_phylum <- summary_phylum %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

## by treatment and phylum
summary_treatment_p <- sum_stats(p.df, 
                                 "Treatment", 
                                 "Phylum", 
                                 variable = "total")
# sort
summary_treatment_p <- summary_treatment_p %>%
  arrange(Treatment, desc(mean))

# round
summary_treatment_p <- summary_treatment_p %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

## by family
summary_family <- sum_stats(f.df, 
                            "Family", 
                            variable = "total")
# sort
summary_family <- summary_family %>%
  arrange(desc(mean))

# round
summary_family <- summary_family %>%
  mutate(across(where(is.numeric), ~ round(., 3)))


## by treatment and family
summary_treatment_f <- sum_stats(f.df, 
                                 "Treatment", 
                                 "Family", 
                                 variable = "total")

# sort
summary_treatment_f <- summary_treatment_f %>%
  arrange(Treatment, desc(mean))

# round
summary_treatment_f <- summary_treatment_f %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

## save ##
write.csv(x = summary_phylum, 
          file = "mag_phylum_relative_abundance_summary.csv", 
          row.names = F)

write.csv(x = summary_treatment_p, 
          file = "mag_phylum_treatment_relative_abundance_summary.csv", 
          row.names = F)

write.csv(x = summary_family, 
          file = "mag_family_relative_abundance_summary.csv", 
          row.names = F)

write.csv(x = summary_treatment_f, 
          file = "mag_family_treatment_relative_abundance_summary.csv", 
          row.names = F)












