# Project:
# Purpose: Fit vital rate models to test the effect of grass-endophyte symbiosis and endophyte hyphal density on  vital rate models (survival, growth, flowering and spikelet).
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

# names(demography_climate)
# view(demography_climate)
# summary(demography_climate)
# Diagnostic of the response variable
ggplot(demography_climate_distance, aes(x = spikelet_t1, fill = Species, color = Species)) +
  geom_density(alpha = 0.4) + # Density plot with transparency
  labs(
    x = "Spikelet average", y = "Density",
    title = ""
  ) +
  theme_bw() +
  # scale_color_manual(values = c("red", "#00AFBB")) +  # Custom colors for the species
  # scale_fill_manual(values = c("red", "#00AFBB")) +
  theme(
    legend.position = c(0.8, 0.8), # Legend at the top
    legend.background = element_rect(fill = "white", color = "white"), # Optional: outline the legend
    legend.title = element_text(size = 10), # Optional: adjust legend title size
    legend.text = element_text(size = 8) # Optional: adjust legend text size
  )

ggplot(demography_climate_distance, aes(x = grow, fill = Species, color = Species)) +
  geom_density(alpha = 0.4) + # Density plot with transparency
  labs(
    x = "Relative growth", y = "Density",
    title = ""
  ) +
  theme_bw() +
  # scale_color_manual(values = c("red", "#00AFBB")) +  # Custom colors for the species
  # scale_fill_manual(values = c("red", "#00AFBB")) +
  theme(
    legend.position = c(0.8, 0.8), # Legend at the top
    legend.background = element_rect(fill = "white", color = "white"), # Optional: outline the legend
    legend.title = element_text(size = 10), # Optional: adjust legend title size
    legend.text = element_text(size = 8) # Optional: adjust legend text size
  )
# Plot the data with regression lines
filtered_data <- demography_climate_distance %>% filter(Species %in% c("ELVI", "POAU"))
ggplot(filtered_data, aes(x = spikelet_t1, y = inf_t1, color = Species)) +
  geom_point() + # Plot points
  geom_smooth(method = "lm", se = TRUE, aes(color = Species)) + # Add regression line
  labs(x = "Spikelet average", y = "Inflorescence average", title = "") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))

# Probability density plot for Inflorescence grouped by species
ggplot(filtered_data, aes(x = inf_t1, fill = Species, color = Species)) +
  geom_density(alpha = 0.4) + # Density plot with transparency
  labs(
    x = "Inflorescence average", y = "Density",
    title = ""
  ) +
  theme_bw() +
  # scale_color_manual(values = c("red", "#00AFBB")) +  # Custom colors for the species
  # scale_fill_manual(values = c("red", "#00AFBB")) +
  theme(
    legend.position = c(0.8, 0.8), # Legend at the top
    legend.background = element_rect(fill = "white", color = "black"), # Optional: outline the legend
    legend.title = element_text(size = 10), # Optional: adjust legend title size
    legend.text = element_text(size = 8) # Optional: adjust legend text size
  )

# I am not sure what prior to use.So I am going to simulate some data and found that the difference
# Define the range of tau values
x <- seq(0.001, 5, length.out = 1000)
# Compute densities for different priors
df <- data.frame(
  x = rep(x, 4),
  density = c(
    dinvgamma(x, 0.1, 0.1),
    dinvgamma(x, 2, 1),
    dhalfcauchy(x, scale = 1),
    dnorm(x, mean = 0, sd = 1) * 2
  ), # Half-Normal (mirrored normal)
  Prior = rep(c(
    "Inv-Gamma(0.1, 0.1)",
    "Inv-Gamma(2, 1)",
    "Half-Cauchy(1)",
    "Half-Normal(1)"
  ), each = length(x))
)

# Plot the priors
ggplot(df, aes(x, density, color = Prior)) +
  geom_line(size = 1) +
  theme_bw() +
  labs(
    title = "",
    x = "tau",
    y = "Density",
    color = "Prior"
  ) +
  theme(legend.position = c(0.8, 0.8))
# The inverse gamma(0.1, 0.1) prior has a very heavy tail, meaning it gives high probability to very large values, which can cause instability.
# The inverse gamma(2, 1) is more reasonable, providing some regularization while still allowing moderate values.
# The half-Cauchy(1) prior has a fat tail but is more controlled compared to inv_gamma(0.1, 0.1).
# The half-Normal(1) prior is much more concentrated near small values, leading to stronger regularization.

## Running the stan model
sim_pars <- list(
  warmup = 1000,
  iter = 4000,
  thin = 2,
  chains = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

# Survival----
## Read and format survival data to build the model
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

## Separate each variable to use the same model stan
### Cumulative precipitation
demography_surv_ppt <- list(
  nSpp = demography_climate_distance_surv$Species %>% n_distinct(),
  nSite = demography_climate_distance_surv$Site %>% n_distinct(),
  nPop = demography_climate_distance_surv$Population %>% n_distinct(),
  nPlot = demography_climate_distance_surv$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_surv$Species,
  site = demography_climate_distance_surv$Site,
  pop= demography_climate_distance_surv$Population,
  plot = demography_climate_distance_surv$site_species_plot,
  clim = as.vector(demography_climate_distance_surv$ppt),
  endo = demography_climate_distance_surv$Endo,
  herb = demography_climate_distance_surv$Herbivory,
  size = demography_climate_distance_surv$log_size_t0,
  y = demography_climate_distance_surv$surv_t1,
  N = nrow(demography_climate_distance_surv)
)

fit_surv_ppt <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/survival.stan",
  data = demography_surv_ppt,
  warmup = sim_pars$warmup,
  seed = 13,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)
summary(fit_surv_ppt)$summary[, c("Rhat", "n_eff")]
posterior_surv_ppt <- as.array(fit_surv_ppt) # Converts to an array
bayesplot::mcmc_trace(posterior_surv_ppt,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3],
                        bclim2[1], bclim2[2], bclim2[3],
                        bendoclim2[1], bendoclim2[2], bendoclim2[3]
                      )
) + theme_bw()
# Create a new data frame for generating predictions
climate_range <- seq(min(demography_surv_ppt$clim), max(demography_surv_ppt$clim), length.out = 30)
endo_status <- c(0, 1) # Endophyte negative and positive
herb_status <- c(0, 1) # Herbivory no and yes
species <- 1:3 # Species 1, 2, 3

# Create a data frame with all combinations
predictions <- expand.grid(clim = climate_range, endo = endo_status, herb = herb_status, species = species)

# Extract posterior samples
posterior_samples <- rstan::extract(fit_surv_ppt)

# Function to calculate predictions based on the posterior samples
get_predictions <- function(clim, endo, herb, species_index, posterior_samples) {
  b0 <- posterior_samples$b0[, species_index]
  bendo <- posterior_samples$bendo[, species_index]
  bherb <- posterior_samples$bherb[, species_index]
  bclim <- posterior_samples$bclim[, species_index]
  bendoclim <- posterior_samples$bendoclim[, species_index]
  bendoherb <- posterior_samples$bendoherb[, species_index]
  bclim2 <- posterior_samples$bclim2[, species_index]
  bendoclim2 <- posterior_samples$bendoclim2[, species_index]
  # Predicted survival (logit scale)
  logit_preds <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb +
    bclim2 * clim^2 +
    bendoclim2 * endo * clim^2
  # Convert logit to probability using logistic function
  pred_probs <- 1 / (1 + exp(-logit_preds))
  return(pred_probs)
}

# Apply the function to generate predictions for all combinations
n_posterior_samples <- length(posterior_samples$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probs_matrix <- matrix(NA, nrow = nrow(predictions), ncol = n_posterior_samples)

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions)) {
  pred_probs_matrix[i, ] <- get_predictions(
    predictions$clim[i], 
    predictions$endo[i], 
    predictions$herb[i], 
    predictions$species[i], 
    posterior_samples
  )
}
species_1_preds <- get_predictions(0.5, 1, 0, 1, posterior_samples)
species_2_preds <- get_predictions(0.2, 0, 1, 2, posterior_samples)
species_3_preds <- get_predictions(-0.3, 1, 1, 3, posterior_samples)

# Convert the matrix into a data frame with the correct structure
pred_probs_df <- as.data.frame(pred_probs_matrix)
colnames(pred_probs_df) <- paste("Posterior_Sample", 1:n_posterior_samples)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probs_df <- cbind(predictions, pred_probs_df)

# Reshape the data frame so we have long format for ggplot
pred_probs_long_df <- gather(pred_probs_df, key = "Posterior_Sample", value = "Pred_Survival", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervals <- pred_probs_long_df %>%
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
observed_data <- data.frame(
  clim = demography_surv_ppt$clim,  # Your climate data
  endo = demography_surv_ppt$endo,  # Your endophyte status data
  herb = demography_surv_ppt$herb,  # Your herbivory status data
  species = demography_surv_ppt$Spp,  # Your species data
  y = demography_surv_ppt$y  # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
ggplot(cred_intervals, aes(x = exp(clim), y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_data, aes(x = exp(clim), y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "Precipitation (mm)",
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

### SPEI (Standardised precipitation-evapotranspiration index)
demography_surv_spei <- list(
  nSpp = demography_climate_distance_surv$Species %>% n_distinct(),
  nSite = demography_climate_distance_surv$Site %>% n_distinct(),
  nPop = demography_climate_distance_surv$Population %>% n_distinct(),
  nPlot = demography_climate_distance_surv$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_surv$Species,
  site = demography_climate_distance_surv$Site,
  pop = demography_climate_distance_surv$Population,
  plot = demography_climate_distance_surv$site_species_plot,
  clim = as.vector(demography_climate_distance_surv$spei),
  endo = demography_climate_distance_surv$Endo,
  herb = demography_climate_distance_surv$Herbivory,
  size = demography_climate_distance_surv$log_size_t0,
  y = demography_climate_distance_surv$surv_t1,
  N = nrow(demography_climate_distance_surv)
)

fit_surv_spei <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/survival.stan",
  data = demography_surv_spei,
  warmup = sim_pars$warmup,
  seed = 13,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)

summary(fit_surv_spei)$summary[, c("Rhat", "n_eff")]
posterior_surv_spei <- as.array(fit_surv_spei) # Converts to an array
bayesplot::mcmc_trace(posterior_surv_spei,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3],
                        bclim2[1], bclim2[2], bclim2[3],
                        bendoclim2[1], bendoclim2[2], bendoclim2[3]
                      )
) + theme_bw()

# Create a data frame with all combinations
predictions_spei <- expand.grid(clim = seq(min(demography_surv_spei$clim), max(demography_surv_spei$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:3 )
# Extract posterior samples
posterior_samples_spei <- rstan::extract(fit_surv_spei)
# Apply the function to generate predictions for all combinations
n_posterior_samples_spei <- length(posterior_samples_spei$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probs_matrix_spei <- matrix(NA, nrow = nrow(predictions_spei), ncol = n_posterior_samples_spei)
# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions_spei)) {
  pred_probs_matrix_spei[i, ] <- get_predictions(
    predictions_spei$clim[i], 
    predictions_spei$endo[i], 
    predictions_spei$herb[i], 
    predictions_spei$species[i], 
    posterior_samples_spei
  )
}

# Convert the matrix into a data frame with the correct structure
pred_probs_spei <- as.data.frame(pred_probs_matrix_spei)
colnames(pred_probs_spei) <- paste("Posterior_Sample", 1:n_posterior_samples_spei)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probs_spei <- cbind(predictions_spei, pred_probs_spei)

# Reshape the data frame so we have long format for ggplot
pred_probs_long_spei <- gather(pred_probs_spei, key = "Posterior_Sample", value = "Pred_Survival", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervals_spei <- pred_probs_long_spei %>%
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
observed_spei <- data.frame(
  clim = demography_surv_spei$clim,  #  climate data
  endo = demography_surv_spei$endo,  #  endophyte status data
  herb = demography_surv_spei$herb,  #  herbivory status data
  species = demography_surv_spei$Spp,  #  species data
  y = demography_surv_spei$y  # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
ggplot(cred_intervals_spei, aes(x = clim, y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_spei, aes(x = clim, y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "SPEI",
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

### Distance from niche centroid
demography_surv_distance <- list(
  nSpp = demography_climate_distance_surv$Species %>% n_distinct(),
  nSite = demography_climate_distance_surv$Site %>% n_distinct(),
  nPop = demography_climate_distance_surv$Population %>% n_distinct(),
  # survival data
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

fit_surv_distance <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/survival_distance.stan",
  data = demography_surv_distance,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)

summary(fit_surv_distance)$summary[, c("Rhat", "n_eff")]
posterior_surv_distance <- as.array(fit_surv_distance) # Converts to an array
bayesplot::mcmc_trace(posterior_surv_distance,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3]
                      )
) + theme_bw()

predictions_distance <- expand.grid(clim = seq(min(demography_surv_distance$clim), max(demography_surv_distance$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:3 )
# Extract posterior samples
posterior_samples_distance <- rstan::extract(fit_surv_distance)
# Apply the function to generate predictions for all combinations
n_posterior_samples_distance <- length(posterior_samples_distance$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probs_matrix_distance <- matrix(NA, nrow = nrow(predictions_distance), ncol = n_posterior_samples_distance)
# Function to calculate predictions based on the posterior samples
get_predictions <- function(clim, endo, herb, species_index, posterior_samples) {
  b0 <- posterior_samples$b0[, species_index]
  bendo <- posterior_samples$bendo[, species_index]
  bherb <- posterior_samples$bherb[, species_index]
  bclim <- posterior_samples$bclim[, species_index]
  bendoclim <- posterior_samples$bendoclim[, species_index]
  bendoherb <- posterior_samples$bendoherb[, species_index]
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
  pred_probs_matrix_distance[i, ] <- get_predictions(
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

## Save RDS file for further use
# saveRDS(fit_surv_ppt, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo model output/fit_surv_ppt.rds')
# saveRDS(fit_surv_spei, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo model output/fit_surv_spei.rds')
# saveRDS(fit_surv_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo model output/fit_surv_distance.rds')
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

## Separate each variable to use the same model stan
### Precipitation
demography_grow_ppt <- list(
  nSpp = demography_climate_distance_grow$Species %>% n_distinct(),
  nSite = demography_climate_distance_grow$Site %>% n_distinct(),
  nPop = demography_climate_distance_grow$Population %>% n_distinct(),
  nPlot = demography_climate_distance_grow$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_grow$Species,
  site = demography_climate_distance_grow$Site,
  pop = demography_climate_distance_grow$Population,
  plot = demography_climate_distance_grow$Plot,
  clim = as.vector(demography_climate_distance_grow$ppt),
  endo = demography_climate_distance_grow$Endo,
  herb = demography_climate_distance_grow$Herbivory,
  size = demography_climate_distance_grow$log_size_t0,
  y = demography_climate_distance_grow$grow,
  N = nrow(demography_climate_distance_grow)
)

fit_grow_ppt <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/growth.stan",
  data = demography_grow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control)

summary(fit_grow_ppt)$summary[, c("Rhat", "n_eff")]
posterior_grow_ppt <- as.array(fit_grow_ppt) # Converts to an array
bayesplot::mcmc_trace(posterior_grow_ppt,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3],
                        bclim2[1], bclim2[2], bclim2[3],
                        bendoclim2[1], bendoclim2[2], bendoclim2[3]
                      )
) + theme_bw()


# Extract posterior samples
posterior_samples <- rstan::extract(fit_grow_ppt)

# Function to calculate predictions based on the posterior samples
get_predictions_grow <- function(clim, endo, herb, species_index, posterior_samples) {
  b0 <- posterior_samples$b0[, species_index]
  bendo <- posterior_samples$bendo[, species_index]
  bherb <- posterior_samples$bherb[, species_index]
  bclim <- posterior_samples$bclim[, species_index]
  bendoclim <- posterior_samples$bendoclim[, species_index]
  bendoherb <- posterior_samples$bendoherb[, species_index]
  bclim2 <- posterior_samples$bclim2[, species_index]
  bendoclim2 <- posterior_samples$bendoclim2[, species_index]
  # Predicted growth
  predg <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb +
    bclim2 * clim^2 +
    bendoclim2 * endo * clim^2
  # Keep predg
  pred_probs <- predg
  return(pred_probs)
}

# Apply the function to generate predictions for all combinations
n_posterior_samples <- length(posterior_samples$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probg_matrix <- matrix(NA, nrow = nrow(predictions), ncol = n_posterior_samples)

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions)) {
  pred_probg_matrix[i, ] <- get_predictions_grow(
    predictions$clim[i], 
    predictions$endo[i], 
    predictions$herb[i], 
    predictions$species[i], 
    posterior_samples
  )
}
species_1_predg <- get_predictions_grow(0.5, 1, 0, 1, posterior_samples)
species_2_predg <- get_predictions_grow(0.2, 0, 1, 2, posterior_samples)
species_3_predg <- get_predictions_grow(-0.3, 1, 1, 3, posterior_samples)

# Convert the matrix into a data frame with the correct structure
pred_probg_df <- as.data.frame(pred_probg_matrix)
colnames(pred_probg_df) <- paste("Posterior_Sample", 1:n_posterior_samples)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probg_df <- cbind(predictions, pred_probg_df)

# Reshape the data frame so we have long format for ggplot
pred_probg_long_df <- gather(pred_probg_df, key = "Posterior_Sample", value = "Pred_Growth", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervalg <- pred_probg_long_df %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Growth, 0.05),
    upper_90 = quantile(Pred_Growth, 0.95),
    lower_95 = quantile(Pred_Growth, 0.025),
    upper_95 = quantile(Pred_Growth, 0.975),
    median = quantile(Pred_Growth, 0.5),
    mean = mean(Pred_Growth)  # Calculate the mean growth
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed survival)
observed_grow <- data.frame(
  clim = demography_grow_ppt$clim,  # Your climate data
  endo = demography_grow_ppt$endo,  # Your endophyte status data
  herb = demography_grow_ppt$herb,  # Your herbivory status data
  species = demography_grow_ppt$Spp,  # Your species data
  y = demography_grow_ppt$y  # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
ggplot(cred_intervalg, aes(x = exp(clim), y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_grow, aes(x = exp(clim), y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "Precipitation (mm)",
    y = "Predicted realtive growth",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.92, 0.085),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )

### SPEI
demography_grow_spei <- list(
  nSpp = demography_climate_distance_grow$Species %>% n_distinct(),
  nSite = demography_climate_distance_grow$Site %>% n_distinct(),
  nPop = demography_climate_distance_grow$Population %>% n_distinct(),
  nPlot = demography_climate_distance_grow$Plot %>% n_distinct(),
  Spp = demography_climate_distance_grow$Species,
  site = demography_climate_distance_grow$Site,
  pop = demography_climate_distance_grow$Population,
  plot = demography_climate_distance_grow$Plot,
  clim = as.vector(demography_climate_distance_grow$spei),
  endo = demography_climate_distance_grow$Endo,
  herb = demography_climate_distance_grow$Herbivory,
  size = demography_climate_distance_grow$log_size_t0,
  y = demography_climate_distance_grow$grow,
  N = nrow(demography_climate_distance_grow)
)

fit_grow_spei <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/growth.stan",
  data = demography_grow_spei,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control)

summary(fit_grow_spei)$summary[, c("Rhat", "n_eff")]
posterior_grow_spei <- as.array(fit_grow_spei) # Converts to an array
bayesplot::mcmc_trace(posterior_grow_spei,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3],
                        bclim2[1], bclim2[2], bclim2[3],
                        bendoclim2[1], bendoclim2[2], bendoclim2[3]
                      )
) + theme_bw()

# Create a data frame with all combinations
predictions_spei <- expand.grid(clim = seq(min(demography_grow_spei$clim), max(demography_grow_spei$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:3 )
# Extract posterior samples
posterior_samples_spei <- rstan::extract(fit_grow_spei)
# Apply the function to generate predictions for all combinations
n_posterior_samples_spei <- length(posterior_samples_spei$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probg_matrix_spei <- matrix(NA, nrow = nrow(predictions_spei), ncol = n_posterior_samples_spei)
# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions_spei)) {
  pred_probg_matrix_spei[i, ] <- get_predictions(
    predictions_spei$clim[i], 
    predictions_spei$endo[i], 
    predictions_spei$herb[i], 
    predictions_spei$species[i], 
    posterior_samples_spei
  )
}

# Convert the matrix into a data frame with the correct structure
pred_probg_spei <- as.data.frame(pred_probg_matrix_spei)
colnames(pred_probg_spei) <- paste("Posterior_Sample", 1:n_posterior_samples_spei)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probg_spei <- cbind(predictions_spei, pred_probg_spei)

# Reshape the data frame so we have long format for ggplot
pred_probg_long_spei <- gather(pred_probg_spei, key = "Posterior_Sample", value = "Pred_Growth", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean growival probability
cred_intervalg_spei <- pred_probg_long_spei %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Growth, 0.05),
    upper_90 = quantile(Pred_Growth, 0.95),
    lower_95 = quantile(Pred_Growth, 0.025),
    upper_95 = quantile(Pred_Growth, 0.975),
    median = quantile(Pred_Growth, 0.5),
    mean = mean(Pred_Growth)  # Calculate the mean growival probability
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed growival)
observed_spei <- data.frame(
  clim = demography_grow_spei$clim,  #  climate data
  endo = demography_grow_spei$endo,  #  endophyte status data
  herb = demography_grow_spei$herb,  #  herbivory status data
  species = demography_grow_spei$Spp,  #  species data
  y = demography_grow_spei$y  # Observed growival
)

# Plot the results with credible intervals, mean growival, and observed points using ggplot2
ggplot(cred_intervalg_spei, aes(x = clim, y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median growival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean growival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_spei, aes(x = clim, y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "SPEI",
    y = "Predicted relative growth",
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

### Distance fro niche centroid
demography_grow_distance <- list(
  nSpp= demography_climate_distance_grow$Species %>% n_distinct(),
  nSite = demography_climate_distance_grow$Site %>% n_distinct(),
  nPop = demography_climate_distance_grow$Population %>% n_distinct(),
  nPlot = demography_climate_distance_grow$Plot %>% n_distinct(),
  Spp = demography_climate_distance_grow$Species,
  site = demography_climate_distance_grow$Site,
  pop = demography_climate_distance_grow$Population,
  plot = demography_climate_distance_grow$Plot,
  clim = as.vector(demography_climate_distance_grow$distance),
  endo = demography_climate_distance_grow$Endo,
  herb = demography_climate_distance_grow$Herbivory,
  size = demography_climate_distance_grow$log_size_t0,
  y = demography_climate_distance_grow$grow,
  N = nrow(demography_climate_distance_grow)
)

fit_grow_distance <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/growth.stan",
  data = demography_grow_distance,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control)

summary(fit_grow_distance)$summary[, c("Rhat", "n_eff")]
posterior_grow_distance <- as.array(fit_grow_distance) # Converts to an array
bayesplot::mcmc_trace(posterior_grow_distance,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3]
                      )
) + theme_bw()

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
  # Predicted survival (logit scale)
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

# observed_data should have columns: clim, endo, herb, species, y (observed survival)
observed_distance <- data.frame(
  clim = demography_grow_distance$clim,  #  climate data
  endo = demography_grow_distance$endo,  #  endophyte status data
  herb = demography_grow_distance$herb,  #  herbivory status data
  species = demography_grow_distance$Spp,  #  species data
  y = demography_grow_distance$y  # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
ggplot(cred_intervals_distance, aes(x = exp(clim), y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_distance, aes(x = exp(clim), y = y, color = factor(endo)), size = 3) +  # Observed data points
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
    legend.position = c(0.9, 0.8),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
## Save RDS file for further use
# saveRDS(fit_grow_ppt, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_grow_ppt.rds')
# saveRDS(fit_grow_spei, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_grow_spei.rds')
# saveRDS(fit_grow_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_grow_distance.rds')
