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
fit_surv_geo_distance <- readRDS(url("https://www.dropbox.com/scl/fi/ah7h4cn9l9p6gwl3cl0aw/fit_surv_geo_distance.rds?rlkey=dypxz2231pwle7p8mb534uqwr&dl=1"))


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

### Distance from niche center
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

### Distance from geographic center
posterior_samples_surv_geo_distance <- rstan::extract(fit_surv_geo_distance)
# Convert to data frame
posterior_samples_surv_geo_distance_df <- as.data.frame(posterior_samples_surv_geo_distance)
# Get the number of species
n_species <- length(posterior_samples_surv_geo_distance$bendo_s[1, ])
# Convert each coefficient into a long-format data frame
surv_geo_distance_coef_list <- c("b0","bendo", "bherb", "bclim", "bendoclim", "bendoherb")
surv_geo_distance_long_data <- list()
for (coef in surv_geo_distance_coef_list) {
  # Extract the coefficient matrix for the current parameter
  surv_geo_distance_coef_matrix <- posterior_samples_surv_geo_distance[[coef]]
  # Convert to long format
  surv_geo_distance_long_data[[coef]] <- as.data.frame(surv_geo_distance_coef_matrix) %>%
    pivot_longer(cols = everything(), names_to = "species", values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'surv_geo_distance_coef_list'
}
# Combine all into one dataframe
plot_data_surv_geo_distance <- bind_rows(surv_geo_distance_long_data)
# Convert species index to numeric
plot_data_surv_geo_distance$species <- as.numeric(gsub("V", "", plot_data_surv_geo_distance$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_surv_geo_distance <- plot_data_surv_geo_distance %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()
# Change species names
summary_stats_surv_geo_distance$species <- recode(summary_stats_surv_geo_distance$species,
                                              "1" = "AGHY",
                                              "2" = "ELVI",
                                              "3" = "PAOU"
)
# unique(summary_stats_surv_geo_distance$parameter)
# Create the coefficient plot with error bars (credible intervals)
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/geodistance_surv_coeff.pdf", width = 5, height = 9)
ggplot(summary_stats_surv_geo_distance, aes(x = factor(species), y = mean_estimate, color = species)) +
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
  labs(x = "Species", y = "Coefficient estimate", title = "") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("AGHY" = "#00AFBB", "ELVI" = "#E7B800", "PAOU" = "#FC4E07")) # Assign unique colors for species
dev.off()



## Growth----
fit_grow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/85pzzrgvkxwogoybr004t/fit_grow_ppt.rds?rlkey=rz2xlg00u1aqhkxeaix7wiu9e&dl=1"))
fit_grow_spei <- readRDS(url("https://www.dropbox.com/scl/fi/4iuz5ay461qkjv732b5yb/fit_grow_spei.rds?rlkey=gq68ixyetb3v4o3ds7uo24cwk&dl=1"))
fit_grow_distance <- readRDS(url("https://www.dropbox.com/scl/fi/fhwjeizspvvd2dbz11195/fit_grow_distance.rds?rlkey=yt863sjawv1zliwem6a9ved1h&dl=1"))
fit_grow_geo_distance <- readRDS(url("https://www.dropbox.com/scl/fi/m3r1hjc7sv3xvuvue79tz/fit_grow_geo_distance.rds?rlkey=5rj6htaf8fovv2ve8qz0f6hvg&dl=1"))

posterior_samples_grow_ppt <- rstan::extract(fit_grow_ppt)
# Convert to data frame
posterior_samples_grow_ppt_df <- as.data.frame(posterior_samples_grow_ppt)
# Get the number of species
n_species <- length(posterior_samples_grow_ppt$bendo[1, ])
# Convert each coefficient into a long-format data frame
grow_ppt_coef_list <- c("bendo", "bherb", "bclim", "bendoclim", "bendoherb", "bclim2", "bendoclim2")
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
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/geodistance_grow_coeff.pdf", width = 5, height = 9)
ggplot(summary_stats_grow_geo_distance, aes(x = factor(species), y = mean_estimate, color = species)) +
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
  labs(x = "Species", y = "Coefficient estimate", title = "") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("AGHY" = "#00AFBB", "ELVI" = "#E7B800", "PAOU" = "#FC4E07")) # Assign unique colors for species
dev.off()

### Distance from geographic center
posterior_samples_grow_geo_distance <- rstan::extract(fit_grow_geo_distance)
# Convert to data frame
posterior_samples_grow_geo_distance_df <- as.data.frame(posterior_samples_grow_geo_distance)
# Get the number of species
n_species <- length(posterior_samples_grow_geo_distance$bendo[1, ])
# Convert each coefficient into a long-format data frame
grow_geo_distance_coef_list <- c("b0","bendo", "bherb", "bclim", "bendoclim", "bendoherb")
grow_geo_distance_long_data <- list()
for (coef in grow_geo_distance_coef_list) {
  # Extract the coefficient matrix for the current parameter
  grow_geo_distance_coef_matrix <- posterior_samples_grow_geo_distance[[coef]]
  # Convert to long format
  grow_geo_distance_long_data[[coef]] <- as.data.frame(grow_geo_distance_coef_matrix) %>%
    pivot_longer(cols = everything(), names_to = "species", values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'grow_geo_distance_coef_list'
}
# Combine all into one dataframe
plot_data_grow_geo_distance <- bind_rows(grow_geo_distance_long_data)
# Convert species index to numeric
plot_data_grow_geo_distance$species <- as.numeric(gsub("V", "", plot_data_grow_geo_distance$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_grow_geo_distance <- plot_data_grow_geo_distance %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()
# Change species names
summary_stats_grow_geo_distance$species <- recode(summary_stats_grow_geo_distance$species,
                                                  "1" = "AGHY",
                                                  "2" = "ELVI",
                                                  "3" = "PAOU"
)

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/geodistance_grow_coeff.pdf", width = 5, height = 9)
ggplot(summary_stats_grow_geo_distance, aes(x = factor(species), y = mean_estimate, color = species)) +
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
  labs(x = "Species", y = "Coefficient estimate", title = "") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("AGHY" = "#00AFBB", "ELVI" = "#E7B800", "PAOU" = "#FC4E07")) # Assign unique colors for species
dev.off()

## Flowering----
fit_flow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/ajo9euxdtxfmle2g006l3/fit_flow_ppt.rds?rlkey=78cknj1yaqs2m5cyicdzg8up4&dl=1"))
fit_flow_spei <- readRDS(url("https://www.dropbox.com/scl/fi/p4t60pr13u5qtt6hnz8kv/fit_flow_spei.rds?rlkey=5on6zasr6c1qt9wzaqux7lk0d&dl=1"))
fit_flow_distance <- readRDS(url("https://www.dropbox.com/scl/fi/y7qvz4pmy2t2j00gqvwqd/fit_flow_distance.rds?rlkey=qe47fg37dpi4z2lu1ur1lsbkh&dl=1"))
fit_flow_geo_distance <- readRDS(url("https://www.dropbox.com/scl/fi/b2dez0gz22p3m2ke1zs7f/fit_flow_geo_distance.rds?rlkey=9fuwynap20617l1rqffjz62wd&dl=1"))

posterior_samples_flow_ppt <- rstan::extract(fit_flow_ppt)
# Convert to data frame
posterior_samples_flow_ppt_df <- as.data.frame(posterior_samples_flow_ppt)
# Get the number of species
n_species <- length(posterior_samples_flow_ppt$bendo[1, ])
# Convert each coefficient into a long-format data frame
flow_ppt_coef_list <- c("bendo", "bherb", "bclim", "bendoclim", "bendoherb", "bclim2", "bendoclim2")
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

### Distance from geographic center
posterior_samples_flow_geo_distance <- rstan::extract(fit_flow_geo_distance)
# Convert to data frame
posterior_samples_flow_geo_distance_df <- as.data.frame(posterior_samples_flow_geo_distance)
# Get the number of species
n_species <- length(posterior_samples_flow_geo_distance$bendo[1, ])
# Convert each coefficient into a long-format data frame
flow_geo_distance_coef_list <- c("b0","bendo", "bherb", "bclim", "bendoclim", "bendoherb")
flow_geo_distance_long_data <- list()
for (coef in flow_geo_distance_coef_list) {
  # Extract the coefficient matrix for the current parameter
  flow_geo_distance_coef_matrix <- posterior_samples_flow_geo_distance[[coef]]
  # Convert to long format
  flow_geo_distance_long_data[[coef]] <- as.data.frame(flow_geo_distance_coef_matrix) %>%
    pivot_longer(cols = everything(), names_to = "species", values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'flow_geo_distance_coef_list'
}
# Combine all into one dataframe
plot_data_flow_geo_distance <- bind_rows(flow_geo_distance_long_data)
# Convert species index to numeric
plot_data_flow_geo_distance$species <- as.numeric(gsub("V", "", plot_data_flow_geo_distance$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_flow_geo_distance <- plot_data_flow_geo_distance %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()
# Change species names
summary_stats_flow_geo_distance$species <- recode(summary_stats_flow_geo_distance$species,
                                                  "1" = "AGHY",
                                                  "2" = "ELVI",
                                                  "3" = "PAOU"
)


pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/geodistance_flow_coeff.pdf", width = 5, height = 9)
ggplot(summary_stats_flow_geo_distance, aes(x = factor(species), y = mean_estimate, color = species)) +
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
  labs(x = "Species", y = "Coefficient estimate", title = "") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("AGHY" = "#00AFBB", "ELVI" = "#E7B800", "PAOU" = "#FC4E07")) # Assign unique colors for species
dev.off()


## Spikelet----
fit_allsites_spik_aghy_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/8rtir221u5ml997h9usgf/fit_allsites_spik_aghy_ppt.rds?rlkey=4e4ya6tnhnosqqhwu76hfpb2x&dl=1"))
fit_allsites_spik_aghy_spei <- readRDS(url("https://www.dropbox.com/scl/fi/swhi510v2mhnlb62xrcvo/fit_allsites_spik_aghy_spei.rds?rlkey=s9szgdfhokb7jjn7zxix5t93e&dl=1"))
fit_allsites_spik_aghy_distance <- readRDS(url("https://www.dropbox.com/scl/fi/7apjqo8nris1vlgih4enz/fit_allsites_spik_aghy_distance.rds?rlkey=lhcl48ud0uuve6wetkwx4ov14&dl=1"))

posterior_samples_spik_ppt <- rstan::extract(fit_allsites_spik_aghy_ppt)
