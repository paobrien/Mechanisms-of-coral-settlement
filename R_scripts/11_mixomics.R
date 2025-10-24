# Use mixOmics to run sPLS on gene frequency data to identify genes associated settlement

setwd("~/Documents/R/Microbial_inducers")

library(tidyverse)
library(mixOmics)

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

## fortmat data ----

# transpose so samples are row
row.names(ko.df) <- ko.df$ko_id
ko.df <- t(ko.df[,-1])

# make sample ids row in meta
row.names(meta.df) <- meta.df$Sample_ID

# remove blanks and controls 
meta.df <- meta.df %>%
  filter(Treatment != "negative") %>%
  filter(Treatment != "control")

ko.df <- ko.df %>%
  as.data.frame %>%
  dplyr::mutate(Sample_ID = rownames(ko.df)) %>%
  filter(Sample_ID %in% meta.df$Sample_ID) %>%
  dplyr::select(-Sample_ID)

## create dataframes

# dfav 
meta_d <- meta.df %>%
  filter(Coral == "D_favus")

ko_d <- ko.df %>%
  as.data.frame %>%
  dplyr::mutate(Sample_ID = rownames(ko.df)) %>%
  filter(Sample_ID %in% meta_d$Sample_ID) %>%
  dplyr::select(-Sample_ID)

#easp
meta_e <- meta.df %>%
  filter(Coral == "E_aspera") %>%
  filter(Treatment != "D_2M")

ko_e <- ko.df %>%
  as.data.frame %>%
  dplyr::mutate(Sample_ID = rownames(ko.df)) %>%
  filter(Sample_ID %in% meta_e$Sample_ID) %>%
  dplyr::select(-Sample_ID)

# plob
meta_pl <- meta.df %>%
  filter(Coral == "P_lobata")

ko_pl <- ko.df %>%
  as.data.frame %>%
  dplyr::mutate(Sample_ID = rownames(ko.df)) %>%
  filter(Sample_ID %in% meta_pl$Sample_ID) %>%
  dplyr::select(-Sample_ID)

# psin
meta_ps <- meta.df %>%
  filter(Coral == "P_sinensis")

ko_ps <- ko.df %>%
  as.data.frame %>%
  dplyr::mutate(Sample_ID = rownames(ko.df)) %>%
  filter(Sample_ID %in% meta_ps$Sample_ID) %>%
  dplyr::select(-Sample_ID)

# combine to list
# metadata
meta.l <- list(
  "meta_dfav" = meta_d,
  "meta_easp" = meta_e,
  "meta_plob" = meta_pl,
  "meta_psin" = meta_ps
)

# ko dataframe
ko.l <- list(
  "ko_dfav" = ko_d,
  "ko_easp" = ko_e,
  "ko_plob" = ko_pl,
  "ko_psin" = ko_ps
)

# remove ko's no longer in df
# check dims
lapply(ko.l, dim)
lapply(meta.l, function(x) {
  length(x$Sample_ID)
})

# remove ko's
ko.l <- lapply(ko.l, function(df) {
  df <- df[, colSums(df != 0) > 0]
})

# check dims again
lapply(ko.l, dim)

# filter out features with near-zero variance
# freq cut = ratio of the most prevalent value to the second most prevalent
# unique cut = percent of values that are unique
zero_var_features.l <- lapply(ko.l, function(x) {
  nearZeroVar(x, freqCut = 95/5, uniqueCut = 20)
})

# check how many there are with zero/near zero variance
lapply(zero_var_features.l, function(x) {
  length(x$Position)
})

# identify the column index to be removed (position = column)
to_rm.l <- lapply(zero_var_features.l, function(x) {
  x$Position
})

# remove columns from the data frame
for (i in seq_along(zero_var_features.l)) {
  if(length(zero_var_features.l[[i]]$Position) > 0) 
    ko.l[[i]] <- ko.l[[i]][, -to_rm.l[[i]], drop = FALSE]
}
lapply(ko.l, dim)


## PCA ----

# check data first with PCA to detect outliers and general structure of data
tune.pca.l <- lapply(ko.l, function(x) {
  tune.pca(x, ncomp = 10, scale = TRUE) 
})

for (i in seq_along(tune.pca.l)) {
  print(names(tune.pca.l)[i])
  plot(tune.pca.l[[i]])
  readline(prompt = "Press [Enter] to see the next plot...")
}

# save
ggsave("psin_components_pca.svg",  
       device = "svg", 
       width = 12, 
       height = 12, 
       units = "cm", 
       dpi = 300)

# the 'elbow' in the plot can be used to determine how many components to use
# in the PCA, although difficult to visualise above 3

# make a list do determine PCA components
elbow <- list(
  "dfav" = 4,
  "easp" = 3,
  "plob" = 2,
  "psin" = 2
)

# run the PCA and check the components
pca.ko.l <- list()
for (i in seq_along(ko.l)) {
  pca.ko.l[[i]] <- pca(ko.l[[i]], 
                       ncomp = elbow[[i]], 
                       center = TRUE, 
                       scale = TRUE)
  names(pca.ko.l)[i] <- names(ko.l)[i]
}

#check
for(i in seq_along(pca.ko.l)) {
  cat(names(pca.ko.l)[i])
  cat("\n Proportion explained: \n")
  print(pca.ko.l[[i]]$prop_expl_var$X)
  cat("\n Cumulative variance: \n")
  print(pca.ko.l[[i]]$cum.var)
  cat("\n Variance total: \n")
  print(pca.ko.l[[i]]$var.tot)
  cat("\n")
}

## identify the informative variables

# top variables for component 1 only:
lapply(pca.ko.l, function(x) {
  head(selectVar(x, comp = 1)$value)
})

# plot PCA 
# as list
for (i in seq_along(pca.ko.l)) {
  cat(names(pca.ko.l)[i])
  plotIndiv(pca.ko.l[[i]], 
            ind.names = FALSE,
            group = meta.l[[i]]$Settlement_2, 
            #ellipse = TRUE,
            #centroid = TRUE,
            pch = as.factor(meta.l[[i]]$Treatment), 
            legend = TRUE, 
            title = names(pca.ko.l)[i],
            legend.title = 'Treatment', 
            legend.title.pch = 'Treatment')
  readline(prompt = "Press [Enter] to see the next plot...")
}

# for individual
plotIndiv(pca.ko.l$ko_psin, 
          ind.names = FALSE,
          group = as.factor(meta.l$meta_psin$Settlement), 
          #ellipse = TRUE,
          #centroid = TRUE,
          pch = as.factor(meta.l$meta_psin$Treatment), 
          legend = TRUE, 
          title = 'Psin PCA',
          legend.title = 'Settlement', 
          legend.title.pch = 'Treatment')

# save
ggsave("pca_psin.svg", device = "svg", width = 15, height = 12, units = "cm", dpi = 300)

# biplot: plot both samples and variables
biplot(pca.ko.l$ko_plob, 
       #var.names = FALSE,
       group = meta.l$meta_plob$Settlement,
       pch = as.factor(meta.l$meta_plob$Treatment), 
       
       cutoff = 0.935,
       legend.title = 'Settlement',
       legend.title.pch = 'Treatment')

ggsave("biplot_plob_0935.svg", 
       device = "svg",
       width = 25, 
       height = 15, 
       units = "cm", 
       dpi = 300)


## sparse PCA ----

# reduce the number of variables associated with the components

# choose the number of components and variables
# components were evaluated above (elbow)

# set up a grid of keepX values to test (i.e., number of variables to test)
grid.keepX <- c(seq(5, 20, 5), 30, 40, 50, 60)

# optimise the number of variables to keep on each component using repeated cross-validation
tune.spca.result.l <- list()
for (i in seq_along(ko.l)) {
  cat("tuning: ", names(ko.l)[i])
  tune.spca.result.l[[i]] <- tune.spca(ko.l[[i]], 
                                       ncomp = elbow[[i]], 
                                       folds = 5, 
                                       test.keepX = grid.keepX, 
                                       nrepeat = 20) 
  names(tune.spca.result.l)[i] <- names(ko.l)[i]
}

# save tuning 
saveRDS(tune.spca.result.l, "Results/metagenomics/mixOmics/pca/tune_spca_result_l.rds")
tune.spca.result.l <- readRDS("Results/metagenomics/mixOmics/pca/tune_spca_result_l.rds")

# check results
lapply(tune.spca.result.l, function(x) {
  x$choice.keepX
})

# plot results:
# number of variables to select are indicated on the x-axis, the average correlation 
# between predicted and actual components based on cross-validation is calculated 
# and shown on the y-axis for each component. The optimal number of variables to 
# select per component is assessed via one-sided t− tests and is indicated with a diamond.
for (i in seq_along(tune.spca.result.l)) {
  cat(names(tune.spca.result.l)[i])
  plot(tune.spca.result.l[[i]])
  readline(prompt = "Press [Enter] to see the next plot...")
}

# based on the tuning above, perform the final sPCA where the number of variables 
# to select on each component is specified with the argument keepX

# by default center = TRUE, scale = TRUE
keepX.select.l <- lapply(tune.spca.result.l, function(x) {
  x$choice.keepX
})

final.spca.multi.l <- list()
for (i in seq_along(ko.l)) {
  final.spca.multi.l[[i]] <- spca(ko.l[[i]], 
                             ncomp = elbow[[i]], 
                             keepX = keepX.select.l[[i]])
  names(final.spca.multi.l)[i] <- names(ko.l)[i]
}

# save sPCA (can use this to extact variables using selectVar)
saveRDS(final.spca.multi.l, "Results/metagenomics/mixOmics/pca/final_spca_multi_l.rds")

# proportion of explained variance
lapply(final.spca.multi.l, function(x) {
  x$prop_expl_var$X
})

final.spca.multi.l$ko_all$variates

# some explained variance is lost compared to a full PCA, but the aim of this 
# is to identify key genes driving the variation in the data

# examine the sPCA sample plot
for (i in seq_along(final.spca.multi.l)) {
  plotIndiv(final.spca.multi.l[[i]],
            comp = c(1, 2),   # Specify components to plot
            ind.names = FALSE,
            group = meta.l[[i]]$Settlement,
            pch = as.factor(meta.l[[i]]$Treatment),
            title = names(final.spca.multi.l)[i],
            legend = TRUE, 
            legend.title = 'Settlement', 
            legend.title.pch = "Treatment")
  readline(prompt = "Press [Enter] to see the next plot...")
}

ggsave("sPCA_psin.svg", 
       device = "svg",
       width = 15, 
       height = 12, 
       units = "cm", 
       dpi = 300)

# biplot to show selected genes
biplot(final.spca.multi.l$ko_all, 
       group = meta.l$meta_all$Settlement,
       pch = as.factor(meta.l$meta_all$Treatment),
       pch.size = 4, 
       ind.names.size = 0,
       legend =TRUE,
       legend.title = 'Settlement',
       legend.title.pch = 'Treatment')

ggsave("sBiplot_all.svg", 
       device = "svg",
       width = 21, 
       height = 12, 
       units = "cm", 
       dpi = 300)

# extract the variable names and their positive or negative contribution to a 
# given component
select_vars_c1 <- lapply(final.spca.multi.l, function(x) {
   selectVar(x, comp = 1)$value
})

select_vars_c2 <- lapply(final.spca.multi.l, function(x) {
  selectVar(x, comp = 2)$value
})

# save
saveRDS(select_vars_c1, "Results/metagenomics/mixOmics/pca/selected_variables_c1.rds")
saveRDS(select_vars_c2, "Results/metagenomics/mixOmics/pca/selected_variables_c2.rds")

# loading weights can also be visualised, where variables are ranked from the 
# least important (top) to the most important (bottom) 
for (i in seq_along(final.spca.multi.l)) {
  cat(names(final.spca.multi.l)[i])
  plotLoadings(final.spca.multi.l[[i]], comp = 1)
  readline(prompt = "Press [Enter] to see the next plot...")
}

# to save
pdf("loadings_c2_spca_plob.pdf", width = 10, height = 8)

plotLoadings(final.spca.multi.l$ko_plob, comp = 2)

dev.off()


## PLS ----

# Partial Least Squares - explore the relationship between two continuous datasets.
# maximises the covariance between the latent variables, rather than correlation. 
# simultaneously model multiple response variables as well as handle noisy, correlated variables

# Test PLS to module gene abundance with settlement

# sPLS1 regression - explain one Y variable with a combination of selected X variables (ko abundance)

# get y
y.l <- lapply(meta.l, function(x) {
  x$Percent_settled
})

# Define the ‘best’ number of dimensions to explain the data using a PLS1 model with a large number of components. 
# Some of the outputs from the PLS1 object are then retrieved in the perf() function 
# to calculate the Q2 criterion using repeated 10-fold cross-validation
tune.pls1.l <- list()
for (i in seq_along(ko.l)) {
  tune.pls1.l[[i]] <- pls(X = ko.l[[i]], 
                          Y = y.l[[i]], 
                          ncomp = 5, 
                          mode = 'regression')
  names(tune.pls1.l)[i] <- names(ko.l)[i]
}

Q2.pls1.l <- lapply(tune.pls1.l, function(x) {
  perf(x, 
       validation = 'Mfold', 
       folds = 10, 
       nrepeat = 10, 
       progressBar = TRUE)
})

for (i in seq_along(Q2.pls1.l)) {
  plot(Q2.pls1.l[[i]], criterion = 'Q2')
  readline(prompt = "Press [Enter] to see the next plot...")
}

plot(Q2.pls1.l$ko_psin, criterion = 'Q2') # values below the line indicate the threshold where
                                          # adding new dimensions are no longer valuable. in this case 4

# set number of variables to select in X (predictors)

# set a grid of values - thin at the start, but also restricted to a small 
# number of genes for a parsimonious model
list.keepX <- c(5:10, seq(15, 50, 5), seq(60, 100, 10)) 

# edit y - because the below code doesn't like 0's (bug)
y2.l <- lapply(y.l, function(y) {
  y+0.001
}) 

# inspect the keepX grid
tune.spls1.MAE.l <- list()
for (i in seq_along(ko.l)) {
  tune.spls1.MAE.l[[i]] <- tune.spls(ko.l[[i]], 
                                     y2.l[[i]], 
                                     ncomp = 5, 
                                     test.keepX = list.keepX, 
                                     validation = 'Mfold', 
                                     folds = 10,
                                     nrepeat = 10, 
                                     progressBar = TRUE, 
                                     measure = 'MAE')
  names(tune.spls1.MAE.l)[i] <- names(ko.l)[i]
}
saveRDS(tune.spls1.MAE.l, "Results/metagenomics/mixOmics/pls/tune_spls1_MAE_l.rds")

plot(tune.spls1.MAE.l$ko_psin) # look for component (line) and feature number (diamond) with lowest MAE

# optimal number of variables to select in X based on the MAE criterion
choice.ncomp.l <- lapply(tune.spls1.MAE.l, function(x) {
  x$choice.ncomp$ncomp
})

# stop at choice.ncomp
choice.keepX.l <- list()
for (i in seq_along(tune.spls1.MAE.l)) {
  choice.keepX.l[[i]] <- tune.spls1.MAE.l[[i]]$choice.keepX[1:choice.ncomp.l[[i]]]
  names(choice.keepX.l)[i] <- names(tune.spls1.MAE.l)[i]
}
choice.ncomp.l
choice.keepX.l

# final model with tuned parameters
spls1.l <- list()
for (i in seq_along(ko.l)) {
  spls1.l[[i]] <- spls(ko.l[[i]], 
                       y.l[[i]], 
                       ncomp = choice.ncomp.l[[i]],
                       keepX = choice.keepX.l[[i]], 
                       mode = "regression")
  names(spls1.l)[i] <- names(ko.l)[i]
}

# list of genes selected - save to combine with maaslin results
# component 1
selected_vars_c1.l <- lapply(spls1.l, function(x) {
  selectVar(x, comp = 1)$X$name
})
saveRDS(selected_vars_c1.l, "Results/metagenomics/mixOmics/pls/spls1_selected_variables_c1.rds")

# component 2
selected_vars_c2.l <- lapply(spls1.l, function(x) {
  if (x$ncomp > 1) {
    selectVar(x, comp = 2)$X$name
  }
})
saveRDS(selected_vars_c2.l, "Results/metagenomics/mixOmics/pls/spls1_selected_variables_c2.rds")

# component 3
selected_vars_c3 <- selectVar(spls1.l$ko_plob, comp = 3)$X$name
saveRDS(selected_vars_c3, "Results/metagenomics/mixOmics/pls/spls1_selected_variables_c3.rds")

# component 4
selected_vars_c4 <- selectVar(spls1.l$ko_plob, comp = 4)$X$name
saveRDS(selected_vars_c4, "Results/metagenomics/mixOmics/pls/spls1_selected_variables_c4.rds")

# compare the amount of explained variance for the X data set based on the 
# sPLS1 (on 1 component) versus PLS1
for (i in seq_along(spls1.l)) {
  print(names(spls1.l)[i])
  print(spls1.l[[i]]$prop_expl_var$X)
  print(tune.pls1.l[[i]]$prop_expl_var$X)
}

# plotting - add a second dimension to the model, which can include 
# the same number of keepX variables as in the first dimension. However, the 
# interpretation should primarily focus on the first dimension 
spls1.c2.l <- list()
for (i in seq_along(ko.l)) {
  spls1.c2.l[[i]] <- spls(ko.l[[i]], 
                          y2.l[[i]], 
                          ncomp = 2, 
                          keepX = c(rep(choice.keepX.l[[i]], 2)), # update to fix keepX
                          mode = "regression")
  names(spls1.c2.l)[i] <- names(ko.l)[i]
}

for (i in seq_along(meta.l)) {
  meta.l[[i]]$Percent_settled <- replace(meta.l[[i]]$Percent_settled,
                                         meta.l[[i]]$Percent_settled == 100, 
                                         99.99)
}

for (i in seq_along(spls1.c2.l)) {
  print(names(spls1.c2.l)[i])
  plotIndiv(spls1.c2.l[[i]],
            group = meta.l[[i]]$Percent_settled,
            pch = as.factor(meta.l[[i]]$Treatment),
            legend = TRUE, 
            legend.title = 'Settlement', 
            legend.title.pch = 'Treatment')
  readline(prompt = "Press [Enter] for next plot...")
}

ggsave("spls_psin_all.svg", 
       device = "svg", 
       width = 25, 
       height = 16, 
       units = "cm", 
       dpi = 300)


# X variable (KOs) are in agreement with the Y variable (settlement)
# for component 1. Low settlement on left, high on right 

# alternatively we can plot the component associated to the X data set (here 
# corresponding to a linear combination of the selected genes) vs. the 
# component associated to the y variable


## performance assessment of sPLS1

# performance of the final model can be assessed with the perf() function, using repeated cross-validation (CV)

# PLS1 model and performance
pls1.l <- list()
for (i in seq_along(ko.l)) {
  pls1.l[[i]] <- pls(ko.l[[i]], 
                y.l[[i]], 
                ncomp = choice.ncomp.l[[i]], 
                mode = "regression")
  names(pls1.l)[i] <- names(ko.l)[i]
}

perf.pls1.l <- list()
for (i in seq_along(pls1.l)) {
  print(names(pls1.l)[i])
  perf.pls1.l[[i]] <- perf(pls1.l[[i]], 
                           validation = "Mfold", 
                           folds =10, 
                           nrepeat = 50, 
                           progressBar = TRUE)
  names(perf.pls1.l)[i] <- names(pls1.l)[i]
}

lapply(perf.pls1.l, function(x) {
  x$measures$MSEP$summary
  x$measures$MSEP$val
})
perf.pls1.l$ko_psin$measures$MSEP$values # to see values from all reps


# A lower MSEP (Mean Squared Error of Prediction ) indicates better predictive accuracy 
# interpreted here as the X variables selected (listed above with selectVar()) 
# can be considered as a good linear combination of predictors to explain y

# plot the final gene selection as a heatmap
# plot displays the similarity values (pearson correlation coefficient) 
# between the X and Y variables selected across two dimensions, 
# and clustered with a complete Euclidean distance method.

# set colours for color.mixo (needs numeric vector)
# settlement category
plot_clrs_set <- function(df) {
  df %>%
    mutate(tmp = case_when(
      Settlement_2 == "Low" ~ 1,
      Settlement_2 == "Med" ~ 2,
      Settlement_2 == "High" ~ 3,
      TRUE ~ NA_real_
    )) %>%
    pull(tmp)
}

# treatment
plot_clrs_treat <- function(df) {
  df %>%
    mutate(tmp = case_when(
      Treatment == "L_1M" ~ 1,
      Treatment == "L_2M" ~ 2,
      Treatment == "D_2M" ~ 3,
      TRUE ~ NA_real_
    )) %>%
    pull(tmp)
}

for (i in seq_along(meta.l)) {
  meta.l[[i]]$plot_clrs_set <- plot_clrs_set(meta.l[[i]])
  meta.l[[i]]$plot_clrs_treat <- plot_clrs_treat(meta.l[[i]])
}

# for all components selected in model
# settlement clrs
for (i in seq_along(spls1.l)) {
  print(names(spls1.l[i]))
  cim(spls1.l[[i]], 
      xlab = "Genes", 
      ylab = "Sample",
      row.sideColors = color.mixo(meta.l[[i]]$plot_clrs_set),
      save = 'pdf', 
      name.save = paste0("cim_settlement_", names(spls1.l)[i]))
  readline(prompt = "Press [Enter] to view next plot...")
}

# treatment clrs
for (i in seq_along(spls1.l)) {
  print(names(spls1.l[i]))
  cim(spls1.l[[i]], 
      xlab = "Genes", 
      ylab = "Sample",
      row.sideColors = color.mixo(meta.l[[i]]$plot_clrs_treat),
      save = 'pdf', 
      name.save = paste0("cim_treatment_", names(spls1.l)[i]))
  readline(prompt = "Press [Enter] to view next plot...")
}




