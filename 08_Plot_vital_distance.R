# Project:
# Purpose: Plot  vital rate models (survival, growth, flowering and spikelet) vs distance from niche center.
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
# Define some basic functions that we'll use later
quote_bare <- function(...) {
  substitute(alist(...)) %>%
    eval() %>%
    sapply(deparse)
}
set.seed(13)
# Demographic data -----
# Merge the demographic census
datini <- read.csv("https://www.dropbox.com/scl/fi/exwmw8z8vp1qkf8inyeoq/Initialdata.csv?rlkey=kez08s92dgh9v0i08269kx1iq&dl=1", stringsAsFactors = F)
dat23 <- read.csv("https://www.dropbox.com/scl/fi/9ob0vpu2xdq8x7u48866s/census2023.csv?rlkey=i2loj3fezymq1p41bo5lsj3uj&dl=1", stringsAsFactors = F)
dat24 <- read.csv("https://www.dropbox.com/scl/fi/s8pnf1j7c85g6jwc944vw/census2024.csv?rlkey=kwt2x8k16q4w7gndm42komj6o&dl=1", stringsAsFactors = F)
datherbivory <- read.csv("https://www.dropbox.com/scl/fi/suy4twdhy36el0k7ytqsi/herbivory.csv?rlkey=hs4xbjn1zrpnpitry30ng538d&dl=1", stringsAsFactors = F)
# unique(datini$Site)
# unique(datini$dat23)
# unique(datini$dat24)
# names(dat23)
# calculate the average spikelet and inflorescence number for each census
dat23 %>%
  mutate(spikelet_23 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = T)), digit = 0) -> dat23_spike
dat24 %>%
  mutate(spikelet_24 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = T), digit = 0), Inf_24 = round(rowMeans(across(attachedInf_24:brokenInf_24), na.rm = T), digit = 0)) -> dat24_spike

## Merge the initial data with the 23 data and the 23 data with the 24 -----
datini23 <- left_join(x = datini, y = dat23_spike, by = c("Tag_ID"))
names(datini23)
dat2324 <- left_join(x = datini23, y = dat24_spike, by = c("Tag_ID"))
names(dat2324)
dat2324 %>%
  mutate(
    tiller_t = Tiller_23,
    tiller_t1 = Tiller_24,
    inf_t = Inf_23,
    inf_t1 = Inf_24,
    spikelet_t = spikelet_23,
    spikelet_t1 = spikelet_24,
    tiller_Herb_t = tiller_Herb,
    tiller_Herb_t1 = tiller_herb_24
  ) %>%
  dplyr::select(
    Site,
    Species,
    Plot,
    Position,
    Tag_ID,
    Population,
    Clone,
    GreenhouseID,
    Endo,
    tiller_t,
    tiller_t1,
    inf_t,
    inf_t1,
    spikelet_t,
    spikelet_t1,
    tiller_Herb_t,
    tiller_Herb_t1,
    date_23,
    date_24
  ) -> dat2324_t_t1

## Merge the demographic data with the herbivory data -----
dat2324_t_t1_herb <- left_join(x = dat2324_t_t1, y = datherbivory, by = c("Site", "Plot", "Species")) # Merge the demographic data with the herbivory data
head(dat2324_t_t1_herb)
unique(dat2324_t_t1_herb$Species)
# view(dat2324_t_t1_herb)

dat2324_t_t1_herb %>%
  filter(tiller_t1 > 0) %>%
  dplyr::select(Species, spikelet_t1) %>%
  group_by(Species) %>%
  summarise(n = sum(spikelet_t1, na.rm = T))

# Climatic data ----
climate_summary <- readRDS(url("https://www.dropbox.com/scl/fi/z7a57xv1ago4erqrnp0tx/prism_means.rds?rlkey=z0ddxpr7ls4k0x527k5pp2wsx&dl=1"))
climate_summary %>%
  rename(Site = site) -> climate_site
distance_species <- readRDS(url("https://www.dropbox.com/scl/fi/kv9j0n2pbiqgrfnm5a4wn/distance_species.rds?rlkey=vni9e8tjw9enwki0mwgnllzjc&dl=1"))
distance_species %>%
  rename(Site = site_code) -> distance_species_clean

## Merge the demographic data with the climatic data -----
demography_climate <- left_join(x = dat2324_t_t1_herb, y = climate_site, by = c("Site"))
demography_climate_distance <- left_join(x = demography_climate, y = distance_species_clean, by = c("Site", "Species"))

## Create new variables
demography_climate_distance %>%
  mutate(
    surv1 = 1 * (!is.na(demography_climate_distance$tiller_t) & !is.na(demography_climate_distance$tiller_t1)),
    site_species_plot = interaction(demography_climate_distance$Site, demography_climate_distance$Species, demography_climate_distance$Plot),
    grow = (log(demography_climate_distance$tiller_t1 + 1) - log(demography_climate_distance$tiller_t + 1))
  ) -> demography_climate_distance

# Survival----
## Read and format survival data
demography_climate_distance %>%
  subset(tiller_t > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, surv1, sum_ppt, mean_pet, mean_spei, distance
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo = as.integer(factor(Endo)) - 1,
    Herbivory = as.integer(factor(Herbivory)) - 1
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    surv_t1 = surv1,
    ppt = log(sum_ppt),
    pet = log(mean_pet),
    spei = mean_spei,
    distance = log(distance)
  ) -> demography_climate_distance_surv

### Separate each variable to use the same model stan
### Distance from niche centroid
demography_surv_distance <- list(
  nSpp = demography_climate_distance_surv$Species %>% n_distinct(),
  nSite = demography_climate_distance_surv$Site %>% n_distinct(),
  nPop = demography_climate_distance_surv$Population %>% n_distinct(),
  nPlot = demography_climate_distance_surv$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_surv$Species,
  site = demography_climate_distance_surv$Site,
  pop = demography_climate_distance_surv$Population,
  plot = demography_climate_distance_surv$Plot,
  clim = as.vector(demography_climate_distance_surv$distance),
  endo = demography_climate_distance_surv$Endo,
  herb = demography_climate_distance_surv$Herbivory,
  size = demography_climate_distance_surv$log_size_t0,
  y = demography_climate_distance_surv$surv_t1,
  N = nrow(demography_climate_distance_surv)
)

fit_surv_distance<-readRDS(url("https://www.dropbox.com/scl/fi/gxc8edjzdjvsrtlb8zm5o/fit_surv_distance.rds?rlkey=bmtq9q0bxf6ooafq8pz3tyavp&dl=1"))

predictions_distance <- expand.grid(clim = seq(min(demography_surv_distance$clim), max(demography_surv_distance$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:3 )
# Extract posterior samples
posterior_samples_distance <- rstan::extract(fit_surv_distance)
# Apply the function to generate predictions for all combinations
n_posterior_samples_distance <- length(posterior_samples_distance$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probs_matrix_distance <- matrix(NA, nrow = nrow(predictions_distance), ncol = n_posterior_samples_distance)
# Function to calculate predictions based on the posterior samples
get_predictions_distance <- function(clim, endo, herb, species_index, posterior_samples_distance) {
  b0 <- posterior_samples_distance$b0[, species_index]
  bendo <- posterior_samples_distance$bendo[, species_index]
  bherb <- posterior_samples_distance$bherb[, species_index]
  bclim <- posterior_samples_distance$bclim[, species_index]
  bendoclim <- posterior_samples_distance$bendoclim[, species_index]
  bendoherb <- posterior_samples_distance$bendoherb[, species_index]
  # Predicted survival (logit scale)
  logit_preds <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb 
  # Convert logit to probability using logistic function
  pred_probs <- 1 / (1 + exp(-logit_preds))
  return(pred_probs)
}

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions_distance)) {
  pred_probs_matrix_distance[i, ] <- get_predictions_distance(
    predictions_distance$clim[i], 
    predictions_distance$endo[i], 
    predictions_distance$herb[i], 
    predictions_distance$species[i], 
    posterior_samples_distance
  )
}

# Convert the matrix into a data frame with the correct structure
pred_probs_distance <- as.data.frame(pred_probs_matrix_distance)
colnames(pred_probs_distance) <- paste("Posterior_Sample", 1:n_posterior_samples_distance)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probs_distance <- cbind(predictions_distance, pred_probs_distance)

# Reshape the data frame so we have long format for ggplot
pred_probs_long_distance <- gather(pred_probs_distance, key = "Posterior_Sample", value = "Pred_Survival", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervals_distance <- pred_probs_long_distance %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Survival, 0.05),
    upper_90 = quantile(Pred_Survival, 0.95),
    lower_95 = quantile(Pred_Survival, 0.025),
    upper_95 = quantile(Pred_Survival, 0.975),
    median = quantile(Pred_Survival, 0.5),
    mean = mean(Pred_Survival)  # Calculate the mean survival probability
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed survival)
observed_distance <- data.frame(
  clim = demography_surv_distance$clim,  #  climate data
  endo = demography_surv_distance$endo,  #  endophyte status data
  herb = demography_surv_distance$herb,  #  herbivory status data
  species = demography_surv_distance$Spp,  #  species data
  y = demography_surv_distance$y  # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PrSurival_distance.pdf", useDingbats = F, height = 9, width = 7)
ggplot(cred_intervals_distance, aes(x = exp(clim), y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_distance, aes(x = exp(clim), y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "Mahalanobis distance",
    y = "Predicted survival probability",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.8),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

# Growth----
## Read and format survival data to build the model
demography_climate_distance %>%
  subset(tiller_t > 0 & tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, grow, sum_ppt, mean_pet, mean_spei, distance
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo = as.integer(factor(Endo)) - 1,
    Herbivory = as.integer(factor(Herbivory)) - 1
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    grow = grow,
    ppt = log(sum_ppt),
    pet = log(mean_pet),
    spei = mean_spei,
    distance = log(distance)
  ) -> demography_climate_distance_grow

table(demography_climate_distance_grow$Species)
table(demography_climate_distance_grow$Site)
### Distance fro niche centroid
demography_grow_distance <- list(
  nSpp= demography_climate_distance_grow$Species %>% n_distinct(),
  nSite = demography_climate_distance_grow$Site %>% n_distinct(),
  nPop = demography_climate_distance_grow$Population %>% n_distinct(),
  nPlot = demography_climate_distance_grow$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_grow$Species,
  site = demography_climate_distance_grow$Site,
  pop = demography_climate_distance_grow$Population,
  plot = demography_climate_distance_grow$site_species_plot,
  clim = as.vector(demography_climate_distance_grow$distance),
  endo = demography_climate_distance_grow$Endo,
  herb = demography_climate_distance_grow$Herbivory,
  size = demography_climate_distance_grow$log_size_t0,
  y = demography_climate_distance_grow$grow,
  N = nrow(demography_climate_distance_grow)
)

fit_grow_distance <- readRDS(url("https://www.dropbox.com/scl/fi/fhwjeizspvvd2dbz11195/fit_grow_distance.rds?rlkey=yt863sjawv1zliwem6a9ved1h&dl=1"))

predictions_distance <- expand.grid(clim = seq(min(demography_grow_distance$clim), max(demography_grow_distance$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:3 )
# Extract posterior samples
posterior_samples_distance <- rstan::extract(fit_grow_distance)
# Apply the function to generate predictions for all combinations
n_posterior_samples_distance <- length(posterior_samples_distance$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probg_matrix_distance <- matrix(NA, nrow = nrow(predictions_distance), ncol = n_posterior_samples_distance)
# Function to calculate predictions based on the posterior samples
get_predictions_distance <- function(clim, endo, herb, species_index, posterior_samples_distance) {
  b0 <- posterior_samples_distance$b0[, species_index]
  bendo <- posterior_samples_distance$bendo[, species_index]
  bherb <- posterior_samples_distance$bherb[, species_index]
  bclim <- posterior_samples_distance$bclim[, species_index]
  bendoclim <- posterior_samples_distance$bendoclim[, species_index]
  bendoherb <- posterior_samples_distance$bendoherb[, species_index]
  # Predicted survival
  predg <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb 
  # Convert logit to probability using logistic function
  pred_probg <-predg
  return(pred_probg)
}

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions_distance)) {
  pred_probg_matrix_distance[i, ] <- get_predictions_distance(
    predictions_distance$clim[i], 
    predictions_distance$endo[i], 
    predictions_distance$herb[i], 
    predictions_distance$species[i], 
    posterior_samples_distance
  )
}

# Convert the matrix into a data frame with the correct structure
pred_probg_distance <- as.data.frame(pred_probg_matrix_distance)
colnames(pred_probg_distance) <- paste("Posterior_Sample", 1:n_posterior_samples_distance)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probg_distance <- cbind(predictions_distance, pred_probg_distance)

# Reshape the data frame so we have long format for ggplot
pred_probg_long_distance <- gather(pred_probg_distance, key = "Posterior_Sample", value = "Pred_Growth", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean growth
cred_intervals_distance <- pred_probg_long_distance %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Growth, 0.05),
    upper_90 = quantile(Pred_Growth, 0.95),
    lower_95 = quantile(Pred_Growth, 0.025),
    upper_95 = quantile(Pred_Growth, 0.975),
    median = quantile(Pred_Growth, 0.5),
    mean = mean(Pred_Growth)  # Calculate the mean survival probability
  ) %>%
  ungroup()

sum(is.na(cred_intervals_distance))
sum(is.na(observed_distance))
# observed_data should have columns: clim, endo, herb, species, y (observed survival)
observed_distance <- data.frame(
  clim = demography_grow_distance$clim,  #  climate data
  endo = demography_grow_distance$endo,  #  endophyte status data
  herb = demography_grow_distance$herb,  #  herbivory status data
  species = demography_grow_distance$Spp,  #  species data
  y = demography_grow_distance$y  # Observed survival
)

str(cred_intervals_distance)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
ggplot(cred_intervals_distance, aes(x = clim, y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_distance, aes(x = clim, y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "Mahalanobis distance",
    y = "Predicted relative growth",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.5),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )

unique(cred_intervals_distance$endo)
unique(observed_distance$endo)
