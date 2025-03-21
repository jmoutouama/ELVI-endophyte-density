# Project:
# Purpose: Fit vital rate models to test the effect of grass-endophyte symbiosis and endophyte hyphal density on  vital rate models (survival, growth, flowering,fertility).
# Note: Raster files are too large to provide in public repository. They are stored on a local machine
# Authors: Jacob Moutouama
# Date last modified (Y-M-D):
rm(list = ls())
# load packages
# remove.packages(c("StanHeaders", "rstan"))
# install.packages("rstan", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
library(rstan)
# set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(13)
# Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
options(tidyverse.quiet = TRUE)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
library(bayesplot)
# install.packages("countreg",repos = "http://R-Forge.R-project.org")
# library(countreg)
library(rmutil)
library(actuar)
# library(SPEI)
library(LaplacesDemon)
library(ggpubr)
library(raster)
# library(rgdal)
library(readxl)
library(ggsci)
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("scater")
# library(scater)
library(BiocManager)
library(swfscMisc)
library(bayesplot)
library(extraDistr)

# Surivival models
fit_surv_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/hi11gxhpqlrdfg389ir0w/fit_surv_ppt.rds?rlkey=22ujyjnm74c6pw9biw50uasvu&dl=1"))
fit_surv_spei <- readRDS(url("https://www.dropbox.com/scl/fi/0js0md2myjvl2scu69bnm/fit_surv_spei.rds?rlkey=scn11z3a3epfgis8y8jrxke91&dl=1"))
fit_surv_distance <- readRDS(url("https://www.dropbox.com/scl/fi/gxc8edjzdjvsrtlb8zm5o/fit_surv_distance.rds?rlkey=bmtq9q0bxf6ooafq8pz3tyavp&dl=1"))

## Plot the coefficients
### Precipitation
posterior_samples_surv_ppt <- rstan::extract(fit_surv_ppt)
# Convert to data frame
posterior_samples_surv_ppt_df <- as.data.frame(posterior_samples_surv_ppt)
# Get the number of species
n_species <- length(posterior_samples_surv_ppt$bendo_s[1, ])
# Convert each coefficient into a long-format data frame
surv_ppt_coef_list <- c("b0","bendo", "bherb", "bclim", "bendoclim", "bendoherb", "bclim2", "bendoclim2")
surv_ppt_long_data <- list()
for (coef in surv_ppt_coef_list) {
  # Extract the coefficient matrix for the current parameter
  surv_ppt_coef_matrix <- posterior_samples_surv_ppt[[coef]]
  # Convert to long format
  surv_ppt_long_data[[coef]] <- as.data.frame(surv_ppt_coef_matrix) %>%
    pivot_longer(cols = everything(), names_to = "species", values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'surv_ppt_coef_list'
}
# Combine all into one dataframe
plot_data_surv_ppt <- bind_rows(surv_ppt_long_data)
# Convert species index to numeric
plot_data_surv_ppt$species <- as.numeric(gsub("V", "", plot_data_surv_ppt$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_surv_ppt <- plot_data_surv_ppt %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()
# Change species names
summary_stats_surv_ppt$species <- recode(summary_stats_surv_ppt$species,
                                         "1" = "AGHY",
                                         "2" = "ELVI",
                                         "3" = "PAOU"
)
# unique(summary_stats_surv_ppt$parameter)
# Create the coefficient plot with error bars (credible intervals)
ggplot(summary_stats_surv_ppt, aes(x = factor(species), y = mean_estimate, color = species)) +
  geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI),
                  position = position_dodge(width = 0.6),
                  size = 1
  ) + # Adds the error bars with credible intervals
  facet_grid(parameter ~ ., 
             scales = "free_y",
             labeller = labeller(parameter = as_labeller(
               c("b0" = "Intercept", 
                 "bendo" = "Endophyte",
                 "bherb" = "Herbivory",
                 "bclim" = "Climate",
                 "bendoclim" = "Endophyte:Climate",
                 "bendoherb" = "Endophyte:Herbivory",
                 "bclim2" = "b[clim]^2",  # Use plotmath syntax
                 "bendoclim2" = "b[endo:clim]^2"  # Use plotmath syntax
               ), 
               default = label_parsed  # This tells ggplot to interpret as expressions
             ))) + # Facet by parameter
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Add horizontal dashed line at y = 0
  theme_bw() +
  labs(x = "Species", y = "Coefficient Estimate", title = "") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("AGHY" = "blue", "ELVI" = "red", "PAOU" = "green")) # Assign unique colors for species

### SPEI
posterior_samples_surv_spei <- rstan::extract(fit_surv_spei)
# Convert to data frame
posterior_samples_surv_spei_df <- as.data.frame(posterior_samples_surv_spei)
# Get the number of species
n_species <- length(posterior_samples_surv_spei$bendo_s[1, ])
# Convert each coefficient into a long-format data frame
surv_spei_coef_list <- c("b0","bendo", "bherb", "bclim", "bendoclim", "bendoherb", "bclim2", "bendoclim2")
surv_spei_long_data <- list()
for (coef in surv_spei_coef_list) {
  # Extract the coefficient matrix for the current parameter
  surv_spei_coef_matrix <- posterior_samples_surv_spei[[coef]]
  # Convert to long format
  surv_spei_long_data[[coef]] <- as.data.frame(surv_spei_coef_matrix) %>%
    pivot_longer(cols = everything(), names_to = "species", values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'surv_spei_coef_list'
}
# Combine all into one dataframe
plot_data_surv_spei <- bind_rows(surv_spei_long_data)
# Convert species index to numeric
plot_data_surv_spei$species <- as.numeric(gsub("V", "", plot_data_surv_spei$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_surv_spei <- plot_data_surv_spei %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()
# Change species names
summary_stats_surv_spei$species <- recode(summary_stats_surv_spei$species,
                                         "1" = "AGHY",
                                         "2" = "ELVI",
                                         "3" = "PAOU"
)
# unique(summary_stats_surv_spei$parameter)
# Create the coefficient plot with error bars (credible intervals)
ggplot(summary_stats_surv_spei, aes(x = factor(species), y = mean_estimate, color = species)) +
  geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI),
                  position = position_dodge(width = 0.6),
                  size = 1
  ) + # Adds the error bars with credible intervals
  facet_grid(parameter ~ ., 
             scales = "free_y",
             labeller = labeller(parameter = as_labeller(
               c("b0" = "Intercept", 
                 "bendo" = "Endophyte",
                 "bherb" = "Herbivory",
                 "bclim" = "Climate",
                 "bendoclim" = "Endophyte:Climate",
                 "bendoherb" = "Endophyte:Herbivory",
                 "bclim2" = "b[clim]^2",  # Use plotmath syntax
                 "bendoclim2" = "b[endo:clim]^2"  # Use plotmath syntax
               ), 
               default = label_parsed  # This tells ggplot to interpret as expressions
             ))) + # Facet by parameter
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Add horizontal dashed line at y = 0
  theme_bw() +
  labs(x = "Species", y = "Coefficient Estimate", title = "") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("AGHY" = "blue", "ELVI" = "red", "PAOU" = "green")) # Assign unique colors for species

### Distance
posterior_samples_surv_distance <- rstan::extract(fit_surv_distance)
# Convert to data frame
posterior_samples_surv_distance_df <- as.data.frame(posterior_samples_surv_distance)
# Get the number of species
n_species <- length(posterior_samples_surv_distance$bendo_s[1, ])
# Convert each coefficient into a long-format data frame
surv_distance_coef_list <- c("b0","bendo", "bherb", "bclim", "bendoclim", "bendoherb")
surv_distance_long_data <- list()
for (coef in surv_distance_coef_list) {
  # Extract the coefficient matrix for the current parameter
  surv_distance_coef_matrix <- posterior_samples_surv_distance[[coef]]
  # Convert to long format
  surv_distance_long_data[[coef]] <- as.data.frame(surv_distance_coef_matrix) %>%
    pivot_longer(cols = everything(), names_to = "species", values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'surv_distance_coef_list'
}
# Combine all into one dataframe
plot_data_surv_distance <- bind_rows(surv_distance_long_data)
# Convert species index to numeric
plot_data_surv_distance$species <- as.numeric(gsub("V", "", plot_data_surv_distance$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_surv_distance <- plot_data_surv_distance %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()
# Change species names
summary_stats_surv_distance$species <- recode(summary_stats_surv_distance$species,
                                          "1" = "AGHY",
                                          "2" = "ELVI",
                                          "3" = "PAOU"
)
# unique(summary_stats_surv_distance$parameter)
# Create the coefficient plot with error bars (credible intervals)
ggplot(summary_stats_surv_distance, aes(x = factor(species), y = mean_estimate, color = species)) +
  geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI),
                  position = position_dodge(width = 0.6),
                  size = 1
  ) + # Adds the error bars with credible intervals
  facet_grid(parameter ~ ., 
             scales = "free_y",
             labeller = labeller(parameter = as_labeller(
               c("b0" = "Intercept", 
                 "bendo" = "Endophyte",
                 "bherb" = "Herbivory",
                 "bclim" = "Climate",
                 "bendoclim" = "Endophyte:Climate",
                 "bendoherb" = "Endophyte:Herbivory"
               ), 
               default = label_parsed  # This tells ggplot to interpret as expressions
             ))) + # Facet by parameter
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Add horizontal dashed line at y = 0
  theme_bw() +
  labs(x = "Species", y = "Coefficient Estimate", title = "") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("AGHY" = "blue", "ELVI" = "red", "PAOU" = "green")) # Assign unique colors for species




## Growth----
fit_allsites_grow_aghy_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/dqhqkqev0y7uem3jwenc8/fit_allsites_grow_aghy_ppt.rds?rlkey=639eok3pace9i05xaft2dyvbe&dl=1"))
fit_allsites_grow_aghy_pet <- readRDS(url("https://www.dropbox.com/scl/fi/xudvskieehqlyg1i0hljr/fit_allsites_grow_aghy_pet.rds?rlkey=wgk79j044fsumqtb2xpqvwzb1&dl=1"))
fit_allsites_grow_aghy_spei <- readRDS(url("https://www.dropbox.com/scl/fi/ojs1ut85cc650j3o8kub1/fit_allsites_grow_aghy_spei.rds?rlkey=chlnweg1fdql1wshfhwfs7act&dl=1"))
fit_allsites_grow_aghy_distance <- readRDS(url("https://www.dropbox.com/scl/fi/fss05g51j0llees9srzr8/fit_allsites_grow_aghy_distance.rds?rlkey=ni7wzy958fwu10jfsq4svle5r&dl=1"))
fit_allsites_grow_aghy_distance_linear <- readRDS(url("https://www.dropbox.com/scl/fi/cw4dtk21twaillwpb64ul/fit_allsites_grow_aghy_distance_linear.rds?rlkey=chptnx9km9yrmqgikvlvwr73s&dl=1"))

posterior_samples_grow_ppt <- rstan::extract(fit_allsites_grow_aghy_ppt)
# Convert to data frame
posterior_samples_grow_ppt_df <- as.data.frame(posterior_samples_grow_ppt)
# Get the number of species
n_species <- length(posterior_samples_grow_ppt$bendo_g[1, ])
# Convert each coefficient into a long-format data frame
grow_ppt_coef_list <- c("bendo_g", "bherb_g", "bclim_g", "bendoclim_g", "bendoherb_g", "bclim2_g", "bendoclim2_g")
grow_ppt_long_data <- list()
for (coef in grow_ppt_coef_list) {
  # Extract the coefficient matrix for the current parameter
  grow_ppt_coef_matrix <- posterior_samples_grow_ppt[[coef]]
  # Convert to long format
  grow_ppt_long_data[[coef]] <- as.data.frame(grow_ppt_coef_matrix) %>%
    pivot_longer(cols = everything(), names_to = "species", values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'grow_ppt_coef_list'
}

# Combine all into one dataframe
plot_data_grow_ppt <- bind_rows(grow_ppt_long_data)
# Convert species index to numeric
plot_data_grow_ppt$species <- as.numeric(gsub("V", "", plot_data_grow_ppt$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_grow_ppt <- plot_data_grow_ppt %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()

# Change species names
summary_stats_grow_ppt$species <- recode(summary_stats_grow_ppt$species,
                                         "1" = "AGHY",
                                         "2" = "ELVI",
                                         "3" = "PAOU"
)

unique(summary_stats_grow_ppt$parameter)
# Create the coefficient plot with error bars (credible intervals)
ggplot(summary_stats_grow_ppt, aes(x = factor(species), y = mean_estimate, color = species)) +
  geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI),
                  position = position_dodge(width = 0.6),
                  size = 1
  ) + # Adds the error bars with credible intervals
  facet_grid(parameter ~ ., scales = "free") + # Facet by both species and parameter
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Add horizontal dashed line at y = 0
  theme_bw() +
  labs(x = "Species", y = "Coefficient Estimate", title = "") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("AGHY" = "blue", "ELVI" = "red", "PAOU" = "green")) # Assign unique colors for species

## Flowering----
fit_allsites_flow_aghy_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/is4m1b6je38i68dfn6bkf/fit_allsites_flow_aghy_ppt.rds?rlkey=oc3jeswl8jpfkrjb0r20jfnng&dl=1"))
fit_allsites_flow_aghy_pet <- readRDS(url("https://www.dropbox.com/scl/fi/8u4tr13xs7jennhsrnbpq/fit_allsites_flow_aghy_pet.rds?rlkey=f2e65i0lu0j81abyiwwtxl2ss&dl=1"))
fit_allsites_flow_aghy_spei <- readRDS(url("https://www.dropbox.com/scl/fi/yaybwm6qn65vg0swkhing/fit_allsites_flow_aghy_spei.rds?rlkey=ncu9cicj1j658r7usohgvxhcw&dl=1"))
fit_allsites_flow_aghy_distance <- readRDS(url("https://www.dropbox.com/scl/fi/lufyihvzkf7gbl8tfyl1n/fit_allsites_flow_aghy_distance.rds?rlkey=aknvr6uyc5ntyf9m033ycrzrx&dl=1"))
fit_allsites_flow_aghy_distance_linear <- readRDS(url("https://www.dropbox.com/scl/fi/bxcpb1vacfdr2tq368arv/fit_allsites_flow_aghy_distance_linear.rds?rlkey=13v9pfkt59ybfjlyntm8r54ef&dl=1"))

posterior_samples_flow_ppt <- rstan::extract(fit_allsites_flow_aghy_ppt)
# Convert to data frame
posterior_samples_flow_ppt_df <- as.data.frame(posterior_samples_flow_ppt)
# Get the number of species
n_species <- length(posterior_samples_flow_ppt$bendo_f[1, ])
# Convert each coefficient into a long-format data frame
flow_ppt_coef_list <- c("bendo_f", "bherb_f", "bclim_f", "bendoclim_f", "bendoherb_f", "bclim2_f", "bendoclim2_f")
flow_ppt_long_data <- list()
for (coef in flow_ppt_coef_list) {
  # Extract the coefficient matrix for the current parameter
  flow_ppt_coef_matrix <- posterior_samples_flow_ppt[[coef]]
  # Convert to long format
  flow_ppt_long_data[[coef]] <- as.data.frame(flow_ppt_coef_matrix) %>%
    pivot_longer(cols = everything(), names_to = "species", values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'flow_ppt_coef_list'
}

# Combine all into one dataframe
plot_data_flow_ppt <- bind_rows(flow_ppt_long_data)
# Convert species index to numeric
plot_data_flow_ppt$species <- as.numeric(gsub("V", "", plot_data_flow_ppt$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_flow_ppt <- plot_data_flow_ppt %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()

# Change species names
summary_stats_flow_ppt$species <- recode(summary_stats_flow_ppt$species,
                                         "1" = "AGHY",
                                         "2" = "ELVI",
                                         "3" = "PAOU"
)

unique(summary_stats_flow_ppt$parameter)
# Create the coefficient plot with error bars (credible intervals)
ggplot(summary_stats_flow_ppt, aes(x = factor(species), y = mean_estimate, color = species)) +
  geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI),
                  position = position_dodge(width = 0.6),
                  size = 1
  ) + # Adds the error bars with credible intervals
  facet_grid(parameter ~ ., scales = "free") + # Facet by both species and parameter
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Add horizontal dashed line at y = 0
  theme_bw() +
  labs(x = "Species", y = "Coefficient Estimate", title = "") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("AGHY" = "blue", "ELVI" = "red", "PAOU" = "green")) # Assign unique colors for species


## Spikelet----
fit_allsites_spik_aghy_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/8rtir221u5ml997h9usgf/fit_allsites_spik_aghy_ppt.rds?rlkey=4e4ya6tnhnosqqhwu76hfpb2x&dl=1"))
fit_allsites_spik_aghy_pet <- readRDS(url("https://www.dropbox.com/scl/fi/wmn81q56ya2ykf0hk4rrg/fit_allsites_spik_aghy_pet.rds?rlkey=1ifut2cb3zdhh19qm8t53mefh&dl=1"))
fit_allsites_spik_aghy_spei <- readRDS(url("https://www.dropbox.com/scl/fi/swhi510v2mhnlb62xrcvo/fit_allsites_spik_aghy_spei.rds?rlkey=s9szgdfhokb7jjn7zxix5t93e&dl=1"))
fit_allsites_spik_aghy_distance <- readRDS(url("https://www.dropbox.com/scl/fi/7apjqo8nris1vlgih4enz/fit_allsites_spik_aghy_distance.rds?rlkey=lhcl48ud0uuve6wetkwx4ov14&dl=1"))

posterior_samples_spik_ppt <- rstan::extract(fit_allsites_spik_aghy_ppt)
