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
# Define some basic functions that we'll use later
quote_bare <- function(...) {
  substitute(alist(...)) %>%
    eval() %>%
    sapply(deparse)
}

invlogit <- function(x) {
  exp(x) / (1 + exp(x))
}

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots <- length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
      ncol = cols, nrow = ceiling(numPlots / cols)
    )
  }

  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(
        layout.pos.row = matchidx$row,
        layout.pos.col = matchidx$col
      ))
    }
  }
}

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
    site_plot = interaction(demography_climate_distance$Site, demography_climate_distance$Plot),
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

# The inverse gamma(0.1, 0.1) prior has a very heavy tail, meaning it gives high probability to very large values, which can cause instability.
# The inverse gamma(2, 1) is more reasonable, providing some regularization while still allowing moderate values.
# The half-Cauchy(1) prior has a fat tail but is more controlled compared to inv_gamma(0.1, 0.1).
# The half-Normal(1) prior is much more concentrated near small values, leading to stronger regularization.

# Survival----
## Read and format survival data to build the model
demography_climate_distance %>%
  subset(tiller_t > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot, site_plot, Endo, Herbivory,
    tiller_t, surv1, mean_ppt, mean_pet, mean_spei, distance
  ) %>%
  na.omit() %>%
  mutate(
    Site = Site %>% as.factor() %>% as.numeric(),
    Species = Species %>% as.factor() %>% as.numeric(),
    Plot = Plot %>% as.factor() %>% as.numeric(),
    site_plot = site_plot %>% as.factor() %>% as.numeric(),
    Endo = Endo %>% as.factor() %>% as.numeric(),
    Herbivory = Herbivory %>% as.factor() %>% as.numeric(),
    Population = Population %>% as.factor() %>% as.numeric()
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    surv_t1 = surv1,
    ppt = log(mean_ppt),
    pet = log(mean_pet),
    spei = mean_spei,
    distance = log(distance)
  ) -> demography_climate_distance_surv

## Separate each variable to use the same model stan
demography_surv_aghy_ppt <- list(
  n_species = demography_climate_distance_surv$Species %>% n_distinct(),
  n_sites = demography_climate_distance_surv$Site %>% n_distinct(),
  n_pops = demography_climate_distance_surv$Population %>% n_distinct(),
  # survival data
  n_plot_s = demography_climate_distance_surv$Plot %>% n_distinct(),
  species_s = demography_climate_distance_surv$Species,
  site_s = demography_climate_distance_surv$Site,
  pop_s = demography_climate_distance_surv$Population,
  plot_s = demography_climate_distance_surv$Plot,
  clim_s = as.vector(demography_climate_distance_surv$ppt),
  endo_s = demography_climate_distance_surv$Endo - 1,
  herb_s = demography_climate_distance_surv$Herbivory - 1,
  size_s = demography_climate_distance_surv$log_size_t0,
  y_s = demography_climate_distance_surv$surv_t1,
  n_s = nrow(demography_climate_distance_surv)
)

data_sites_surv_aghy_pet <- list(
  n_species = demography_climate_distance_surv$Species %>% n_distinct(),
  n_sites = demography_climate_distance_surv$Site %>% n_distinct(),
  n_pops = demography_climate_distance_surv$Population %>% n_distinct(),
  # survival data
  n_plot_s = demography_climate_distance_surv$Plot %>% n_distinct(),
  species_s = demography_climate_distance_surv$Species,
  site_s = demography_climate_distance_surv$Site,
  pop_s = demography_climate_distance_surv$Population,
  plot_s = demography_climate_distance_surv$Plot,
  clim_s = as.vector(demography_climate_distance_surv$pet),
  endo_s = demography_climate_distance_surv$Endo - 1,
  herb_s = demography_climate_distance_surv$Herbivory - 1,
  size_s = demography_climate_distance_surv$log_size_t0,
  y_s = demography_climate_distance_surv$surv_t1,
  n_s = nrow(demography_climate_distance_surv)
)
data_sites_surv_aghy_spei <- list(
  n_species = demography_climate_distance_surv$Species %>% n_distinct(),
  n_sites = demography_climate_distance_surv$Site %>% n_distinct(),
  n_pops = demography_climate_distance_surv$Population %>% n_distinct(),
  # survival data
  n_plot_s = demography_climate_distance_surv$Plot %>% n_distinct(),
  species_s = demography_climate_distance_surv$Species,
  site_s = demography_climate_distance_surv$Site,
  pop_s = demography_climate_distance_surv$Population,
  plot_s = demography_climate_distance_surv$Plot,
  clim_s = as.vector(demography_climate_distance_surv$spei),
  endo_s = demography_climate_distance_surv$Endo - 1,
  herb_s = demography_climate_distance_surv$Herbivory - 1,
  size_s = demography_climate_distance_surv$log_size_t0,
  y_s = demography_climate_distance_surv$surv_t1,
  n_s = nrow(demography_climate_distance_surv)
)
data_sites_surv_aghy_distance <- list(
  n_species = demography_climate_distance_surv$Species %>% n_distinct(),
  n_sites = demography_climate_distance_surv$Site %>% n_distinct(),
  n_pops = demography_climate_distance_surv$Population %>% n_distinct(),
  # survival data
  n_plot_s = demography_climate_distance_surv$Plot %>% n_distinct(),
  species_s = demography_climate_distance_surv$Species,
  site_s = demography_climate_distance_surv$Site,
  pop_s = demography_climate_distance_surv$Population,
  plot_s = demography_climate_distance_surv$Plot,
  clim_s = as.vector(demography_climate_distance_surv$distance),
  endo_s = demography_climate_distance_surv$Endo - 1,
  herb_s = demography_climate_distance_surv$Herbivory - 1,
  size_s = demography_climate_distance_surv$log_size_t0,
  y_s = demography_climate_distance_surv$surv_t1,
  n_s = nrow(demography_climate_distance_surv)
)
## Running the stan model
sim_pars <- list(
  warmup = 1000,
  iter = 4000,
  thin = 2,
  chains = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

# fit_allsites_surv_aghy_ppt <- stan(
#  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/survival.stan",
#  data = demography_surv_aghy_ppt,
#  warmup = sim_pars$warmup,
#  seed = 13,
#  iter = sim_pars$iter,
#  thin = sim_pars$thin,
#  chains = sim_pars$chains,
#  control = sim_pars$control)

fit_allsites_surv_aghy_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/kpac617znxe2arebi8wqf/fit_allsites_surv_aghy_ppt.rds?rlkey=8okogsfhm0abkeeus810gsgy8&dl=1"))

## Chain mixing and model convergence
summary(fit_allsites_surv_aghy_ppt)$summary[, c("Rhat", "n_eff")]
posterior_surv_aghy_ppt <- as.array(fit_allsites_surv_aghy_ppt) # Converts to an array
bayesplot::mcmc_trace(posterior_surv_aghy_ppt,
  pars = quote_bare(
    b0_s[1], b0_s[2], b0_s[3],
    bendo_s[1], bendo_s[2], bendo_s[3],
    bherb_s[1], bherb_s[2], bherb_s[3],
    bclim_s[1], bclim_s[2], bclim_s[3],
    bendoclim_s[1], bendoclim_s[2], bendoclim_s[3],
    bendoherb_s[1], bendoherb_s[2], bendoherb_s[3],
    bclim2_s[1], bclim2_s[2], bclim2_s[3],
    bendoclim2_s[1], bendoclim2_s[2], bendoclim2_s[3]
  )
) + theme_bw()

## Compute predictions using posterior draws
posterior_samples_sur_ppt <- rstan::extract(fit_allsites_surv_aghy_ppt)
n_draws <- 1000 # Number of posterior samples to use
pred_data_sur_ppt <- as.data.frame(demography_surv_aghy_ppt)
pred_matrix_sur_ppt <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_sur_ppt))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_sur_ppt)) {
    species <- pred_data_sur_ppt$species_s[j]
    pred_matrix_sur_ppt[i, j] <- posterior_samples_sur_ppt$b0_s[i, species] +
      posterior_samples_sur_ppt$bendo_s[i, species] * pred_data_sur_ppt$endo_s[j] +
      posterior_samples_sur_ppt$bclim_s[i, species] * pred_data_sur_ppt$clim_s[j] +
      posterior_samples_sur_ppt$bherb_s[i, species] * pred_data_sur_ppt$herb_s[j] +
      posterior_samples_sur_ppt$bendoclim_s[i, species] * pred_data_sur_ppt$clim_s[j] * pred_data_sur_ppt$endo_s[j] +
      posterior_samples_sur_ppt$bendoherb_s[i, species] * pred_data_sur_ppt$endo_s[j] * pred_data_sur_ppt$herb_s[j] +
      posterior_samples_sur_ppt$bclim2_s[i, species] * pred_data_sur_ppt$clim_s[j]^2 +
      posterior_samples_sur_ppt$bendoclim2_s[i, species] * pred_data_sur_ppt$endo_s[j] * pred_data_sur_ppt$clim_s[j]^2
  }
}

# Convert logits to probability scale
pred_prob_sur_ppt <- invlogit(pred_matrix_sur_ppt)
# Compute mean and 95% credible interval
pred_data_sur_ppt$mean_survival <- apply(pred_prob_sur_ppt, 2, mean)
pred_data_sur_ppt$lower_ci <- apply(pred_prob_sur_ppt, 2, quantile, probs = 0.025)
pred_data_sur_ppt$upper_ci <- apply(pred_prob_sur_ppt, 2, quantile, probs = 0.975)
# Plot predicted survival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PrSurv_ppt.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_sur_ppt, aes(x = exp(clim_s), y = mean_survival, color = factor(endo_s), fill = factor(endo_s))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_s)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_s), y = y_s, color = factor(endo_s)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  facet_grid(species_s ~ herb_s, labeller = labeller(
    species_s = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_s = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Preciptation (mm)",
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

# fit_allsites_surv_aghy_pet <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/survival.stan",
#   data = data_sites_surv_aghy_pet,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control)

fit_allsites_surv_aghy_pet <- readRDS(url("https://www.dropbox.com/scl/fi/e98tc65ecu5kq0z90tmuc/fit_allsites_surv_aghy_pet.rds?rlkey=zav1ijycjc0odvfszc3i17w2a&dl=1"))

summary(fit_allsites_surv_aghy_pet)$summary[, c("Rhat", "n_eff")]
posterior_surv_aghy_pet <- as.array(fit_allsites_surv_aghy_pet) # Converts to an array
bayesplot::mcmc_trace(posterior_surv_aghy_pet,
  pars = quote_bare(
    b0_s[1], b0_s[2], b0_s[3],
    bendo_s[1], bendo_s[2], bendo_s[3],
    bherb_s[1], bherb_s[2], bherb_s[3],
    bclim_s[1], bclim_s[2], bclim_s[3],
    bendoclim_s[1], bendoclim_s[2], bendoclim_s[3],
    bendoherb_s[1], bendoherb_s[2], bendoherb_s[3],
    bclim2_s[1], bclim2_s[2], bclim2_s[3],
    bendoclim2_s[1], bendoclim2_s[2], bendoclim2_s[3]
  )
) + theme_bw()

## Compute predictions using posterior draws
posterior_samples_sur_pet <- rstan::extract(fit_allsites_surv_aghy_pet)
n_draws <- 1000 # Number of posterior samples to use
pred_data_sur_pet <- as.data.frame(data_sites_surv_aghy_pet)
pred_matrix_sur_pet <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_sur_pet))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_sur_pet)) {
    species <- pred_data_sur_pet$species_s[j]
    pred_matrix_sur_pet[i, j] <- posterior_samples_sur_pet$b0_s[i, species] +
      posterior_samples_sur_pet$bendo_s[i, species] * pred_data_sur_pet$endo_s[j] +
      posterior_samples_sur_pet$bclim_s[i, species] * pred_data_sur_pet$clim_s[j] +
      posterior_samples_sur_pet$bherb_s[i, species] * pred_data_sur_pet$herb_s[j] +
      posterior_samples_sur_pet$bendoclim_s[i, species] * pred_data_sur_pet$clim_s[j] * pred_data_sur_pet$endo_s[j] +
      posterior_samples_sur_pet$bendoherb_s[i, species] * pred_data_sur_pet$endo_s[j] * pred_data_sur_pet$herb_s[j] +
      posterior_samples_sur_pet$bclim2_s[i, species] * pred_data_sur_pet$clim_s[j]^2 +
      posterior_samples_sur_pet$bendoclim2_s[i, species] * pred_data_sur_pet$endo_s[j] * pred_data_sur_pet$clim_s[j]^2
  }
}

# Convert logits to probability scale
pred_prob_sur_pet <- invlogit(pred_matrix_sur_pet)
# Compute mean and 95% credible interval
pred_data_sur_pet$mean_survival <- apply(pred_prob_sur_pet, 2, mean)
pred_data_sur_pet$lower_ci <- apply(pred_prob_sur_pet, 2, quantile, probs = 0.025)
pred_data_sur_pet$upper_ci <- apply(pred_prob_sur_pet, 2, quantile, probs = 0.975)
# Plot predicted survival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PrSurv_pet.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_sur_pet, aes(x = exp(clim_s), y = mean_survival, color = factor(endo_s), fill = factor(endo_s))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_s)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_s), y = y_s, color = factor(endo_s)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  facet_grid(species_s ~ herb_s, labeller = labeller(
    species_s = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_s = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Potential Evapotranspiration",
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

# fit_allsites_surv_aghy_spei <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/survival.stan",
#   data = data_sites_surv_aghy_spei,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control)

fit_allsites_surv_aghy_spei <- readRDS(url("https://www.dropbox.com/scl/fi/0wx8mr8ignd36dl5utlrj/fit_allsites_surv_aghy_spei.rds?rlkey=i1w37cl91823m7eioi2gup80i&dl=1"))

summary(fit_allsites_surv_aghy_spei)$summary[, c("Rhat", "n_eff")]
posterior_surv_aghy_spei <- as.array(fit_allsites_surv_aghy_spei) # Converts to an array
bayesplot::mcmc_trace(posterior_surv_aghy_spei,
  pars = quote_bare(
    b0_s[1], b0_s[2], b0_s[3],
    bendo_s[1], bendo_s[2], bendo_s[3],
    bherb_s[1], bherb_s[2], bherb_s[3],
    bclim_s[1], bclim_s[2], bclim_s[3],
    bendoclim_s[1], bendoclim_s[2], bendoclim_s[3],
    bendoherb_s[1], bendoherb_s[2], bendoherb_s[3],
    bclim2_s[1], bclim2_s[2], bclim2_s[3],
    bendoclim2_s[1], bendoclim2_s[2], bendoclim2_s[3]
  )
) + theme_bw()

## Compute predictions using posterior draws
posterior_samples_sur_spei <- rstan::extract(fit_allsites_surv_aghy_spei)
n_draws <- 1000 # Number of posterior samples to use
pred_data_sur_spei <- as.data.frame(data_sites_surv_aghy_spei)
pred_matrix_sur_spei <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_sur_spei))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_sur_spei)) {
    species <- pred_data_sur_spei$species_s[j]
    pred_matrix_sur_spei[i, j] <- posterior_samples_sur_spei$b0_s[i, species] +
      posterior_samples_sur_spei$bendo_s[i, species] * pred_data_sur_spei$endo_s[j] +
      posterior_samples_sur_spei$bclim_s[i, species] * pred_data_sur_spei$clim_s[j] +
      posterior_samples_sur_spei$bherb_s[i, species] * pred_data_sur_spei$herb_s[j] +
      posterior_samples_sur_spei$bendoclim_s[i, species] * pred_data_sur_spei$clim_s[j] * pred_data_sur_spei$endo_s[j] +
      posterior_samples_sur_spei$bendoherb_s[i, species] * pred_data_sur_spei$endo_s[j] * pred_data_sur_spei$herb_s[j] +
      posterior_samples_sur_spei$bclim2_s[i, species] * pred_data_sur_spei$clim_s[j]^2 +
      posterior_samples_sur_spei$bendoclim2_s[i, species] * pred_data_sur_spei$endo_s[j] * pred_data_sur_spei$clim_s[j]^2
  }
}

# Convert logits to probability scale
pred_prob_sur_spei <- invlogit(pred_matrix_sur_spei)
# Compute mean and 95% credible interval
pred_data_sur_spei$mean_survival <- apply(pred_prob_sur_spei, 2, mean)
pred_data_sur_spei$lower_ci <- apply(pred_prob_sur_spei, 2, quantile, probs = 0.025)
pred_data_sur_spei$upper_ci <- apply(pred_prob_sur_spei, 2, quantile, probs = 0.975)
# Plot predicted survival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PrSurv_spei.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_sur_spei, aes(x = exp(clim_s), y = mean_survival, color = factor(endo_s), fill = factor(endo_s))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_s)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_s), y = y_s, color = factor(endo_s)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  facet_grid(species_s ~ herb_s, labeller = labeller(
    species_s = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_s = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Standardised precipitation-evapotranspiration index",
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
# fit_allsites_surv_aghy_distance <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/survival.stan",
#   data = data_sites_surv_aghy_distance,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control)

fit_allsites_surv_aghy_distance <- readRDS(url("https://www.dropbox.com/scl/fi/ap0uys7p36wy6fl5zm4gk/fit_allsites_surv_aghy_distance.rds?rlkey=ro1rdaplg8hx9dfbx9v5jw605&dl=1"))

summary(fit_allsites_surv_aghy_distance)$summary[, c("Rhat", "n_eff")]
posterior_surv_aghy_distance <- as.array(fit_allsites_surv_aghy_distance) # Converts to an array
bayesplot::mcmc_trace(posterior_surv_aghy_distance,
  pars = quote_bare(
    b0_s[1], b0_s[2], b0_s[3],
    bendo_s[1], bendo_s[2], bendo_s[3],
    bherb_s[1], bherb_s[2], bherb_s[3],
    bclim_s[1], bclim_s[2], bclim_s[3],
    bendoclim_s[1], bendoclim_s[2], bendoclim_s[3],
    bendoherb_s[1], bendoherb_s[2], bendoherb_s[3],
    bclim2_s[1], bclim2_s[2], bclim2_s[3],
    bendoclim2_s[1], bendoclim2_s[2], bendoclim2_s[3]
  )
) + theme_bw()

## Compute predictions using posterior draws
posterior_samples_sur_distance <- rstan::extract(fit_allsites_surv_aghy_distance)
n_draws <- 1000 # Number of posterior samples to use
pred_data_sur_distance <- as.data.frame(data_sites_surv_aghy_distance)
pred_matrix_sur_distance <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_sur_distance))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_sur_distance)) {
    species <- pred_data_sur_distance$species_s[j]
    pred_matrix_sur_distance[i, j] <- posterior_samples_sur_distance$b0_s[i, species] +
      posterior_samples_sur_distance$bendo_s[i, species] * pred_data_sur_distance$endo_s[j] +
      posterior_samples_sur_distance$bclim_s[i, species] * pred_data_sur_distance$clim_s[j] +
      posterior_samples_sur_distance$bherb_s[i, species] * pred_data_sur_distance$herb_s[j] +
      posterior_samples_sur_distance$bendoclim_s[i, species] * pred_data_sur_distance$clim_s[j] * pred_data_sur_distance$endo_s[j] +
      posterior_samples_sur_distance$bendoherb_s[i, species] * pred_data_sur_distance$endo_s[j] * pred_data_sur_distance$herb_s[j] +
      posterior_samples_sur_distance$bclim2_s[i, species] * pred_data_sur_distance$clim_s[j]^2 +
      posterior_samples_sur_distance$bendoclim2_s[i, species] * pred_data_sur_distance$endo_s[j] * pred_data_sur_distance$clim_s[j]^2
  }
}

# Convert logits to probability scale
pred_prob_sur_distance <- invlogit(pred_matrix_sur_distance)
# Compute mean and 95% credible interval
pred_data_sur_distance$mean_survival <- apply(pred_prob_sur_distance, 2, mean)
pred_data_sur_distance$lower_ci <- apply(pred_prob_sur_distance, 2, quantile, probs = 0.025)
pred_data_sur_distance$upper_ci <- apply(pred_prob_sur_distance, 2, quantile, probs = 0.975)
# Plot predicted survival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PrSurv_mh.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_sur_distance, aes(x = exp(clim_s), y = mean_survival, color = factor(endo_s), fill = factor(endo_s))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +  # Credible interval
  # geom_point(aes(x = exp(clim_s), y = y_s, color = factor(endo_s)),
  #            position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5) +
  facet_grid(species_s ~ herb_s, labeller = labeller(
    species_s = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_s = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Mahalanobis distance",
    y = "Predicted survival probability",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "red"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "red"), labels = c("E-", "E+")) + # Change fill labels
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

# fit_allsites_surv_aghy_distance_linear <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/survival_distance.stan",
#   data = data_sites_surv_aghy_distance,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control
# )

fit_allsites_surv_aghy_distance_linear <- readRDS(url("https://www.dropbox.com/scl/fi/3u88o5n2bq1zb74wvsgx7/fit_allsites_surv_aghy_distance_linear.rds?rlkey=igquh02xfwrh1dp8danzv6w0r&dl=1"))

summary(fit_allsites_surv_aghy_distance_linear)$summary[, c("Rhat", "n_eff")]
posterior_surv_aghy_distance_linear <- as.array(fit_allsites_surv_aghy_distance_linear) # Converts to an array
bayesplot::mcmc_trace(posterior_surv_aghy_distance_linear,
  pars = quote_bare(
    b0_s[1], b0_s[2], b0_s[3],
    bendo_s[1], bendo_s[2], bendo_s[3],
    bherb_s[1], bherb_s[2], bherb_s[3],
    bclim_s[1], bclim_s[2], bclim_s[3],
    bendoclim_s[1], bendoclim_s[2], bendoclim_s[3],
    bendoherb_s[1], bendoherb_s[2], bendoherb_s[3]
  )
) + theme_bw()

## Compute predictions using posterior draws
posterior_samples_sur_distance <- rstan::extract(fit_allsites_surv_aghy_distance_linear)
n_draws <- 1000 # Number of posterior samples to use
pred_data_sur_distance <- as.data.frame(data_sites_surv_aghy_distance)
pred_matrix_sur_distance <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_sur_distance))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_sur_distance)) {
    species <- pred_data_sur_distance$species_s[j]
    pred_matrix_sur_distance[i, j] <- posterior_samples_sur_distance$b0_s[i, species] +
      posterior_samples_sur_distance$bendo_s[i, species] * pred_data_sur_distance$endo_s[j] +
      posterior_samples_sur_distance$bclim_s[i, species] * pred_data_sur_distance$clim_s[j] +
      posterior_samples_sur_distance$bherb_s[i, species] * pred_data_sur_distance$herb_s[j] +
      posterior_samples_sur_distance$bendoclim_s[i, species] * pred_data_sur_distance$clim_s[j] * pred_data_sur_distance$endo_s[j] +
      posterior_samples_sur_distance$bendoherb_s[i, species] * pred_data_sur_distance$endo_s[j] * pred_data_sur_distance$herb_s[j]
  }
}

# Convert logits to probability scale
pred_prob_sur_distance <- invlogit(pred_matrix_sur_distance)
# Compute mean and 95% credible interval
pred_data_sur_distance$mean_survival <- apply(pred_prob_sur_distance, 2, mean)
pred_data_sur_distance$lower_ci <- apply(pred_prob_sur_distance, 2, quantile, probs = 0.025)
pred_data_sur_distance$upper_ci <- apply(pred_prob_sur_distance, 2, quantile, probs = 0.975)
# Plot predicted survival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PrSurv_mh_l.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_sur_distance, aes(x = exp(clim_s), y = mean_survival, color = factor(endo_s), fill = factor(endo_s))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +  # Credible interval
  # geom_point(aes(x = exp(clim_s), y = y_s, color = factor(endo_s)),
  #            position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5) +
  facet_grid(species_s ~ herb_s, labeller = labeller(
    species_s = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_s = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Mahalanobis distance",
    y = "Predicted survival probability",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "red"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "red"), labels = c("E-", "E+")) + # Change fill labels
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





## Save RDS file for further use
# saveRDS(fit_allsites_surv_aghy_ppt, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_surv_aghy_ppt.rds')
# saveRDS(fit_allsites_surv_aghy_pet, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_surv_aghy_pet.rds')
# saveRDS(fit_allsites_surv_aghy_spei, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_surv_aghy_spei.rds')
# saveRDS(fit_allsites_surv_aghy_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_surv_aghy_distance.rds')
# saveRDS(fit_allsites_surv_aghy_distance_linear, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_surv_aghy_distance_linear.rds')

# Growth----
## Read and format survival data to build the model
demography_climate_distance %>%
  subset(tiller_t > 0 & tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot, site_plot, Endo, Herbivory,
    tiller_t, grow, mean_ppt, mean_pet, mean_spei, distance
  ) %>%
  na.omit() %>%
  mutate(
    Site = Site %>% as.factor() %>% as.numeric(),
    Species = Species %>% as.factor() %>% as.numeric(),
    Plot = Plot %>% as.factor() %>% as.numeric(),
    site_plot = site_plot %>% as.factor() %>% as.numeric(),
    Endo = Endo %>% as.factor() %>% as.numeric(),
    Herbivory = Herbivory %>% as.factor() %>% as.numeric(),
    Population = Population %>% as.factor() %>% as.numeric()
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    grow = grow,
    ppt = log(mean_ppt),
    pet = log(mean_pet),
    spei = mean_spei,
    distance = log(distance)
  ) -> demography_climate_distance_grow

## Separate each variable to use the same model stan
demography_grow_aghy_ppt <- list(
  n_species = demography_climate_distance_grow$Species %>% n_distinct(),
  n_sites = demography_climate_distance_grow$Site %>% n_distinct(),
  n_pops = demography_climate_distance_grow$Population %>% n_distinct(),
  # growth data
  n_plot_g = demography_climate_distance_grow$Plot %>% n_distinct(),
  species_g = demography_climate_distance_grow$Species,
  site_g = demography_climate_distance_grow$Site,
  pop_g = demography_climate_distance_grow$Population,
  plot_g = demography_climate_distance_grow$Plot,
  clim_g = as.vector(demography_climate_distance_grow$ppt),
  endo_g = demography_climate_distance_grow$Endo - 1,
  herb_g = demography_climate_distance_grow$Herbivory - 1,
  size_g = demography_climate_distance_grow$log_size_t0,
  y_g = demography_climate_distance_grow$grow,
  n_g = nrow(demography_climate_distance_grow)
)

data_sites_grow_aghy_pet <- list(
  n_species = demography_climate_distance_grow$Species %>% n_distinct(),
  n_sites = demography_climate_distance_grow$Site %>% n_distinct(),
  n_pops = demography_climate_distance_grow$Population %>% n_distinct(),
  # growth data
  n_plot_g = demography_climate_distance_grow$Plot %>% n_distinct(),
  species_g = demography_climate_distance_grow$Species,
  site_g = demography_climate_distance_grow$Site,
  pop_g = demography_climate_distance_grow$Population,
  plot_g = demography_climate_distance_grow$Plot,
  clim_g = as.vector(demography_climate_distance_grow$pet),
  endo_g = demography_climate_distance_grow$Endo - 1,
  herb_g = demography_climate_distance_grow$Herbivory - 1,
  size_g = demography_climate_distance_grow$log_size_t0,
  y_g = demography_climate_distance_grow$grow,
  n_g = nrow(demography_climate_distance_grow)
)
data_sites_grow_aghy_spei <- list(
  n_species = demography_climate_distance_grow$Species %>% n_distinct(),
  n_sites = demography_climate_distance_grow$Site %>% n_distinct(),
  n_pops = demography_climate_distance_grow$Population %>% n_distinct(),
  # growth data
  n_plot_g = demography_climate_distance_grow$Plot %>% n_distinct(),
  species_g = demography_climate_distance_grow$Species,
  site_g = demography_climate_distance_grow$Site,
  pop_g = demography_climate_distance_grow$Population,
  plot_g = demography_climate_distance_grow$Plot,
  clim_g = as.vector(demography_climate_distance_grow$spei),
  endo_g = demography_climate_distance_grow$Endo - 1,
  herb_g = demography_climate_distance_grow$Herbivory - 1,
  size_g = demography_climate_distance_grow$log_size_t0,
  y_g = demography_climate_distance_grow$grow,
  n_g = nrow(demography_climate_distance_grow)
)
data_sites_grow_aghy_distance <- list(
  n_species = demography_climate_distance_grow$Species %>% n_distinct(),
  n_sites = demography_climate_distance_grow$Site %>% n_distinct(),
  n_pops = demography_climate_distance_grow$Population %>% n_distinct(),
  # growth data
  n_plot_g = demography_climate_distance_grow$Plot %>% n_distinct(),
  species_g = demography_climate_distance_grow$Species,
  site_g = demography_climate_distance_grow$Site,
  pop_g = demography_climate_distance_grow$Population,
  plot_g = demography_climate_distance_grow$Plot,
  clim_g = as.vector(demography_climate_distance_grow$distance),
  endo_g = demography_climate_distance_grow$Endo - 1,
  herb_g = demography_climate_distance_grow$Herbivory - 1,
  size_g = demography_climate_distance_grow$log_size_t0,
  y_g = demography_climate_distance_grow$grow,
  n_g = nrow(demography_climate_distance_grow)
)

# fit_allsites_grow_aghy_ppt <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/growth.stan",
#   data = demography_grow_aghy_ppt,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control)

fit_allsites_grow_aghy_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/dqhqkqev0y7uem3jwenc8/fit_allsites_grow_aghy_ppt.rds?rlkey=639eok3pace9i05xaft2dyvbe&dl=1"))

summary(fit_allsites_grow_aghy_ppt)$summary[, c("Rhat", "n_eff")]
posterior_grow_aghy_ppt <- as.array(fit_allsites_grow_aghy_ppt) # Converts to an array
bayesplot::mcmc_trace(posterior_grow_aghy_ppt,
  pars = quote_bare(
    b0_g[1], b0_g[2], b0_g[3],
    bendo_g[1], bendo_g[2], bendo_g[3],
    bherb_g[1], bherb_g[2], bherb_g[3],
    bclim_g[1], bclim_g[2], bclim_g[3],
    bendoclim_g[1], bendoclim_g[2], bendoclim_g[3],
    bendoherb_g[1], bendoherb_g[2], bendoherb_g[3],
    bclim2_g[1], bclim2_g[2], bclim2_g[3],
    bendoclim2_g[1], bendoclim2_g[2], bendoclim2_g[3]
  )
) + theme_bw()

## Compute predictions using posterior draws
posterior_samples_grow_ppt <- rstan::extract(fit_allsites_grow_aghy_ppt)
n_draws <- 1000 # Number of posterior samples to use
pred_data_grow_ppt <- as.data.frame(demography_grow_aghy_ppt)
pred_matrix_grow_ppt <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_grow_ppt))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_grow_ppt)) {
    species <- pred_data_grow_ppt$species_g[j]
    pred_matrix_grow_ppt[i, j] <- posterior_samples_grow_ppt$b0_g[i, species] +
      posterior_samples_grow_ppt$bendo_g[i, species] * pred_data_grow_ppt$endo_g[j] +
      posterior_samples_grow_ppt$bclim_g[i, species] * pred_data_grow_ppt$clim_g[j] +
      posterior_samples_grow_ppt$bherb_g[i, species] * pred_data_grow_ppt$herb_g[j] +
      posterior_samples_grow_ppt$bendoclim_g[i, species] * pred_data_grow_ppt$clim_g[j] * pred_data_grow_ppt$endo_g[j] +
      posterior_samples_grow_ppt$bendoherb_g[i, species] * pred_data_grow_ppt$endo_g[j] * pred_data_grow_ppt$herb_g[j] +
      posterior_samples_grow_ppt$bclim2_g[i, species] * pred_data_grow_ppt$clim_g[j]^2 +
      posterior_samples_grow_ppt$bendoclim2_g[i, species] * pred_data_grow_ppt$endo_g[j] * pred_data_grow_ppt$clim_g[j]^2
  }
}

# Convert to probability scale
pred_prob_grow_ppt <- exp(pred_matrix_grow_ppt)
# Compute mean and 95% credible interval
pred_data_grow_ppt$mean_growth <- apply(pred_prob_grow_ppt, 2, mean)
pred_data_grow_ppt$lower_ci <- apply(pred_prob_grow_ppt, 2, quantile, probs = 0.025)
pred_data_grow_ppt$upper_ci <- apply(pred_prob_grow_ppt, 2, quantile, probs = 0.975)
# Plot predicted survival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Prgrow_ppt.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_grow_ppt, aes(x = exp(clim_g), y = mean_growth, color = factor(endo_g), fill = factor(endo_g))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_g)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_g), y = y_g, color = factor(endo_g)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  facet_grid(species_g ~ herb_g, labeller = labeller(
    species_g = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_g = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Preciptation (mm)",
    y = "Relative growth",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.90),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()


# fit_allsites_grow_aghy_pet <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/growth.stan",
#   data = data_sites_grow_aghy_pet,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control)

fit_allsites_grow_aghy_pet <- readRDS(url("https://www.dropbox.com/scl/fi/xudvskieehqlyg1i0hljr/fit_allsites_grow_aghy_pet.rds?rlkey=wgk79j044fsumqtb2xpqvwzb1&dl=1"))

summary(fit_allsites_grow_aghy_pet)$summary[, c("Rhat", "n_eff")]
posterior_grow_aghy_pet <- as.array(fit_allsites_grow_aghy_pet) # Converts to an array
bayesplot::mcmc_trace(posterior_grow_aghy_pet,
  pars = quote_bare(
    b0_g[1], b0_g[2], b0_g[3],
    bendo_g[1], bendo_g[2], bendo_g[3],
    bherb_g[1], bherb_g[2], bherb_g[3],
    bclim_g[1], bclim_g[2], bclim_g[3],
    bendoclim_g[1], bendoclim_g[2], bendoclim_g[3],
    bendoherb_g[1], bendoherb_g[2], bendoherb_g[3],
    bclim2_g[1], bclim2_g[2], bclim2_g[3],
    bendoclim2_g[1], bendoclim2_g[2], bendoclim2_g[3]
  )
) + theme_bw()

## Compute predictions using posterior draws
posterior_samples_grow_pet <- rstan::extract(fit_allsites_grow_aghy_pet)
n_draws <- 1000 # Number of posterior samples to use
pred_data_grow_pet <- as.data.frame(data_sites_grow_aghy_pet)
pred_matrix_grow_pet <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_grow_pet))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_grow_pet)) {
    species <- pred_data_grow_pet$species_g[j]
    pred_matrix_grow_pet[i, j] <- posterior_samples_grow_pet$b0_g[i, species] +
      posterior_samples_grow_pet$bendo_g[i, species] * pred_data_grow_pet$endo_g[j] +
      posterior_samples_grow_pet$bclim_g[i, species] * pred_data_grow_pet$clim_g[j] +
      posterior_samples_grow_pet$bherb_g[i, species] * pred_data_grow_pet$herb_g[j] +
      posterior_samples_grow_pet$bendoclim_g[i, species] * pred_data_grow_pet$clim_g[j] * pred_data_grow_pet$endo_g[j] +
      posterior_samples_grow_pet$bendoherb_g[i, species] * pred_data_grow_pet$endo_g[j] * pred_data_grow_pet$herb_g[j] +
      posterior_samples_grow_pet$bclim2_g[i, species] * pred_data_grow_pet$clim_g[j]^2 +
      posterior_samples_grow_pet$bendoclim2_g[i, species] * pred_data_grow_pet$endo_g[j] * pred_data_grow_pet$clim_g[j]^2
  }
}

# Convert logits to probability scale
pred_prob_grow_pet <- exp(pred_matrix_grow_pet)
# Compute mean and 95% credible interval
pred_data_grow_pet$mean_growth <- apply(pred_prob_grow_pet, 2, mean)
pred_data_grow_pet$lower_ci <- apply(pred_prob_grow_pet, 2, quantile, probs = 0.025)
pred_data_grow_pet$upper_ci <- apply(pred_prob_grow_pet, 2, quantile, probs = 0.975)
# Plot predicted growth  with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Prgrow_pet.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_grow_pet, aes(x = exp(clim_g), y = mean_growth, color = factor(endo_g), fill = factor(endo_g))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_g)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_g), y = y_g, color = factor(endo_g)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  facet_grid(species_g ~ herb_g, labeller = labeller(
    species_g = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_g = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Potential Evapotranspiration",
    y = "Predicted relative growth",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.9),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )

dev.off()
# fit_allsites_grow_aghy_spei <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/growth.stan",
#   data = data_sites_grow_aghy_spei,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control)

fit_allsites_grow_aghy_spei <- readRDS(url("https://www.dropbox.com/scl/fi/ojs1ut85cc650j3o8kub1/fit_allsites_grow_aghy_spei.rds?rlkey=chlnweg1fdql1wshfhwfs7act&dl=1"))

summary(fit_allsites_grow_aghy_spei)$summary[, c("Rhat", "n_eff")]
posterior_grow_aghy_spei <- as.array(fit_allsites_grow_aghy_spei) # Converts to an array
bayesplot::mcmc_trace(posterior_grow_aghy_spei,
  pars = quote_bare(
    b0_g[1], b0_g[2], b0_g[3],
    bendo_g[1], bendo_g[2], bendo_g[3],
    bherb_g[1], bherb_g[2], bherb_g[3],
    bclim_g[1], bclim_g[2], bclim_g[3],
    bendoclim_g[1], bendoclim_g[2], bendoclim_g[3],
    bendoherb_g[1], bendoherb_g[2], bendoherb_g[3],
    bclim2_g[1], bclim2_g[2], bclim2_g[3],
    bendoclim2_g[1], bendoclim2_g[2], bendoclim2_g[3]
  )
) + theme_bw()

## Compute predictions using posterior draws
posterior_samples_grow_spei <- rstan::extract(fit_allsites_grow_aghy_spei)
n_draws <- 2000 # Number of posterior samples to use
pred_data_grow_spei <- as.data.frame(data_sites_grow_aghy_spei)
pred_matrix_grow_spei <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_grow_spei))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_grow_spei)) {
    species <- pred_data_grow_spei$species_g[j]
    pred_matrix_grow_spei[i, j] <- posterior_samples_grow_spei$b0_g[i, species] +
      posterior_samples_grow_spei$bendo_g[i, species] * pred_data_grow_spei$endo_g[j] +
      posterior_samples_grow_spei$bclim_g[i, species] * pred_data_grow_spei$clim_g[j] +
      posterior_samples_grow_spei$bherb_g[i, species] * pred_data_grow_spei$herb_g[j] +
      posterior_samples_grow_spei$bendoclim_g[i, species] * pred_data_grow_spei$clim_g[j] * pred_data_grow_spei$endo_g[j] +
      posterior_samples_grow_spei$bendoherb_g[i, species] * pred_data_grow_spei$endo_g[j] * pred_data_grow_spei$herb_g[j] +
      posterior_samples_grow_spei$bclim2_g[i, species] * pred_data_grow_spei$clim_g[j]^2 +
      posterior_samples_grow_spei$bendoclim2_g[i, species] * pred_data_grow_spei$endo_g[j] * pred_data_grow_spei$clim_g[j]^2
  }
}

# Convert logits to probability scale
pred_prob_grow_spei <- exp(pred_matrix_grow_spei)
# Compute mean and 95% credible interval
pred_data_grow_spei$mean_growth <- apply(pred_prob_grow_spei, 2, mean)
pred_data_grow_spei$lower_ci <- apply(pred_prob_grow_spei, 2, quantile, probs = 0.025)
pred_data_grow_spei$upper_ci <- apply(pred_prob_grow_spei, 2, quantile, probs = 0.975)
# Plot predicted growvival probabilities with credible intervals

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Prgrow_spei.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_grow_spei, aes(x = exp(clim_g), y = mean_growth, color = factor(endo_g), fill = factor(endo_g))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_g)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_g), y = y_g, color = factor(endo_g)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  facet_grid(species_g ~ herb_g, labeller = labeller(
    species_g = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_g = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Standardised precipitation-evapotranspiration index",
    y = "Predicted survival probability",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.9),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )

dev.off()
# fit_allsites_grow_aghy_distance <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/growth.stan",
#   data = data_sites_grow_aghy_distance,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control)

fit_allsites_grow_aghy_distance <- readRDS(url("https://www.dropbox.com/scl/fi/fss05g51j0llees9srzr8/fit_allsites_grow_aghy_distance.rds?rlkey=ni7wzy958fwu10jfsq4svle5r&dl=1"))

summary(fit_allsites_grow_aghy_distance)$summary[, c("Rhat", "n_eff")]
posterior_grow_aghy_distance <- as.array(fit_allsites_grow_aghy_distance) # Converts to an array
bayesplot::mcmc_trace(posterior_grow_aghy_distance,
  pars = quote_bare(
    b0_g[1], b0_g[2], b0_g[3],
    bendo_g[1], bendo_g[2], bendo_g[3],
    bherb_g[1], bherb_g[2], bherb_g[3],
    bclim_g[1], bclim_g[2], bclim_g[3],
    bendoclim_g[1], bendoclim_g[2], bendoclim_g[3],
    bendoherb_g[1], bendoherb_g[2], bendoherb_g[3],
    bclim2_g[1], bclim2_g[2], bclim2_g[3],
    bendoclim2_g[1], bendoclim2_g[2], bendoclim2_g[3]
  )
) + theme_bw()


## Compute predictions using posterior draws
posterior_samples_grow_distance <- rstan::extract(fit_allsites_grow_aghy_distance)
n_draws <- 1000 # Number of posterior samples to use
pred_data_grow_distance <- as.data.frame(data_sites_grow_aghy_distance)
pred_matrix_grow_distance <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_grow_distance))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_grow_distance)) {
    species <- pred_data_grow_distance$species_g[j]
    pred_matrix_grow_distance[i, j] <- posterior_samples_grow_distance$b0_g[i, species] +
      posterior_samples_grow_distance$bendo_g[i, species] * pred_data_grow_distance$endo_g[j] +
      posterior_samples_grow_distance$bclim_g[i, species] * pred_data_grow_distance$clim_g[j] +
      posterior_samples_grow_distance$bherb_g[i, species] * pred_data_grow_distance$herb_g[j] +
      posterior_samples_grow_distance$bendoclim_g[i, species] * pred_data_grow_distance$clim_g[j] * pred_data_grow_distance$endo_g[j] +
      posterior_samples_grow_distance$bendoherb_g[i, species] * pred_data_grow_distance$endo_g[j] * pred_data_grow_distance$herb_g[j] +
      posterior_samples_grow_distance$bclim2_g[i, species] * pred_data_grow_distance$clim_g[j]^2 +
      posterior_samples_grow_distance$bendoclim2_g[i, species] * pred_data_grow_distance$endo_g[j] * pred_data_grow_distance$clim_g[j]^2
  }
}

# Convert logits to probability scale
pred_prob_grow_distance <- exp(pred_matrix_grow_distance)
# Compute mean and 95% credible interval
pred_data_grow_distance$mean_growth <- apply(pred_prob_grow_distance, 2, mean)
pred_data_grow_distance$lower_ci <- apply(pred_prob_grow_distance, 2, quantile, probs = 0.025)
pred_data_grow_distance$upper_ci <- apply(pred_prob_grow_distance, 2, quantile, probs = 0.975)
# Plot predicted growvival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Prgrow_mh.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_grow_distance, aes(x = exp(clim_g), y = mean_growth, color = factor(endo_g), fill = factor(endo_g))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_g)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_g), y = y_g, color = factor(endo_g)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  facet_grid(species_g ~ herb_g, labeller = labeller(
    species_g = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_g = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
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
dev.off()

# fit_allsites_grow_aghy_distance_linear <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/growth_distance.stan",
#   data = data_sites_grow_aghy_distance,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control)

fit_allsites_grow_aghy_distance_linear <- readRDS(url("https://www.dropbox.com/scl/fi/cw4dtk21twaillwpb64ul/fit_allsites_grow_aghy_distance_linear.rds?rlkey=chptnx9km9yrmqgikvlvwr73s&dl=1"))

summary(fit_allsites_grow_aghy_distance_linear)$summary[, c("Rhat", "n_eff")]
posterior_grow_aghy_distance <- as.array(fit_allsites_grow_aghy_distance_linear) # Converts to an array
bayesplot::mcmc_trace(posterior_grow_aghy_distance,
  pars = quote_bare(
    b0_g[1], b0_g[2], b0_g[3],
    bendo_g[1], bendo_g[2], bendo_g[3],
    bherb_g[1], bherb_g[2], bherb_g[3],
    bclim_g[1], bclim_g[2], bclim_g[3],
    bendoclim_g[1], bendoclim_g[2], bendoclim_g[3],
    bendoherb_g[1], bendoherb_g[2], bendoherb_g[3]
  )
) + theme_bw()


## Compute predictions using posterior draws
posterior_samples_grow_distance <- rstan::extract(fit_allsites_grow_aghy_distance_linear)
n_draws <- 1000 # Number of posterior samples to use
pred_data_grow_distance <- as.data.frame(data_sites_grow_aghy_distance)
pred_matrix_grow_distance <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_grow_distance))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_grow_distance)) {
    species <- pred_data_grow_distance$species_g[j]
    pred_matrix_grow_distance[i, j] <- posterior_samples_grow_distance$b0_g[i, species] +
      posterior_samples_grow_distance$bendo_g[i, species] * pred_data_grow_distance$endo_g[j] +
      posterior_samples_grow_distance$bclim_g[i, species] * pred_data_grow_distance$clim_g[j] +
      posterior_samples_grow_distance$bherb_g[i, species] * pred_data_grow_distance$herb_g[j] +
      posterior_samples_grow_distance$bendoclim_g[i, species] * pred_data_grow_distance$clim_g[j] * pred_data_grow_distance$endo_g[j] +
      posterior_samples_grow_distance$bendoherb_g[i, species] * pred_data_grow_distance$endo_g[j] * pred_data_grow_distance$herb_g[j]
  }
}

# Convert logits to probability scale
pred_prob_grow_distance <- pred_matrix_grow_distance
# Compute mean and 95% credible interval
pred_data_grow_distance$mean_growth <- apply(pred_prob_grow_distance, 2, mean)
pred_data_grow_distance$lower_ci <- apply(pred_prob_grow_distance, 2, quantile, probs = 0.025)
pred_data_grow_distance$upper_ci <- apply(pred_prob_grow_distance, 2, quantile, probs = 0.975)
# Plot predicted growvival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Prgrow_mh_l.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_grow_distance, aes(x = exp(clim_g), y = mean_growth, color = factor(endo_g), fill = factor(endo_g))) +
  # geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_g)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_g), y = y_g, color = factor(endo_g)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.3) + # Linear regression with confidence interval
  facet_grid(species_g ~ herb_g, labeller = labeller(
    species_g = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_g = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
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
dev.off()



## Save RDS file for further use
# saveRDS(fit_allsites_grow_aghy_ppt, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_grow_aghy_ppt.rds')
# saveRDS(fit_allsites_grow_aghy_pet, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_grow_aghy_pet.rds')
# saveRDS(fit_allsites_grow_aghy_spei, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_grow_aghy_spei.rds')
# saveRDS(fit_allsites_grow_aghy_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_grow_aghy_distance.rds')
# saveRDS(fit_allsites_grow_aghy_distance_linear, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_grow_aghy_distance_linear.rds')

# Flowering----
demography_climate_distance %>%
  subset(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot, site_plot, Endo, Herbivory,
    tiller_t, inf_t1, mean_ppt, mean_pet, mean_spei, distance
  ) %>%
  na.omit() %>%
  mutate(
    Site = Site %>% as.factor() %>% as.numeric(),
    Species = Species %>% as.factor() %>% as.numeric(),
    Plot = Plot %>% as.factor() %>% as.numeric(),
    site_plot = site_plot %>% as.factor() %>% as.numeric(),
    Endo = Endo %>% as.factor() %>% as.numeric(),
    Herbivory = Herbivory %>% as.factor() %>% as.numeric(),
    Population = Population %>% as.factor() %>% as.numeric()
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    flow_t1 = inf_t1,
    ppt = log(mean_ppt),
    pet = log(mean_pet),
    spei = mean_spei,
    distance = log(distance)
  ) -> demography_climate_distance_flow

# hist(demography_climate_distance_flow$flow_t1,main="")
## Separate each variable to use the same model stan
demography_flow_aghy_ppt <- list(
  n_species = demography_climate_distance_flow$Species %>% n_distinct(),
  n_sites = demography_climate_distance_flow$Site %>% n_distinct(),
  n_pops = demography_climate_distance_flow$Population %>% n_distinct(),
  # flowering data
  n_plot_f = demography_climate_distance_flow$Plot %>% n_distinct(),
  species_f = demography_climate_distance_flow$Species,
  site_f = demography_climate_distance_flow$Site,
  pop_f = demography_climate_distance_flow$Population,
  plot_f = demography_climate_distance_flow$Plot,
  clim_f = as.vector(demography_climate_distance_flow$ppt),
  endo_f = demography_climate_distance_flow$Endo - 1,
  herb_f = demography_climate_distance_flow$Herbivory - 1,
  size_f = demography_climate_distance_flow$log_size_t0,
  y_f = demography_climate_distance_flow$flow_t1,
  n_f = nrow(demography_climate_distance_flow)
)

data_sites_flow_aghy_pet <- list(
  n_species = demography_climate_distance_flow$Species %>% n_distinct(),
  n_sites = demography_climate_distance_flow$Site %>% n_distinct(),
  n_pops = demography_climate_distance_flow$Population %>% n_distinct(),
  # flowering data
  n_plot_f = demography_climate_distance_flow$Plot %>% n_distinct(),
  species_f = demography_climate_distance_flow$Species,
  site_f = demography_climate_distance_flow$Site,
  pop_f = demography_climate_distance_flow$Population,
  plot_f = demography_climate_distance_flow$Plot,
  clim_f = as.vector(demography_climate_distance_flow$pet),
  endo_f = demography_climate_distance_flow$Endo - 1,
  herb_f = demography_climate_distance_flow$Herbivory - 1,
  size_f = demography_climate_distance_flow$log_size_t0,
  y_f = demography_climate_distance_flow$flow_t1,
  n_f = nrow(demography_climate_distance_flow)
)
data_sites_flow_aghy_spei <- list(
  n_species = demography_climate_distance_flow$Species %>% n_distinct(),
  n_sites = demography_climate_distance_flow$Site %>% n_distinct(),
  n_pops = demography_climate_distance_flow$Population %>% n_distinct(),
  # flowering data
  n_plot_f = demography_climate_distance_flow$Plot %>% n_distinct(),
  species_f = demography_climate_distance_flow$Species,
  site_f = demography_climate_distance_flow$Site,
  pop_f = demography_climate_distance_flow$Population,
  plot_f = demography_climate_distance_flow$Plot,
  clim_f = as.vector(demography_climate_distance_flow$spei),
  endo_f = demography_climate_distance_flow$Endo - 1,
  herb_f = demography_climate_distance_flow$Herbivory - 1,
  size_f = demography_climate_distance_flow$log_size_t0,
  y_f = demography_climate_distance_flow$flow_t1,
  n_f = nrow(demography_climate_distance_flow)
)
data_sites_flow_aghy_distance <- list(
  n_species = demography_climate_distance_flow$Species %>% n_distinct(),
  n_sites = demography_climate_distance_flow$Site %>% n_distinct(),
  n_pops = demography_climate_distance_flow$Population %>% n_distinct(),
  # flowering data
  n_plot_f = demography_climate_distance_flow$Plot %>% n_distinct(),
  species_f = demography_climate_distance_flow$Species,
  site_f = demography_climate_distance_flow$Site,
  pop_f = demography_climate_distance_flow$Population,
  plot_f = demography_climate_distance_flow$Plot,
  clim_f = as.vector(demography_climate_distance_flow$distance),
  endo_f = demography_climate_distance_flow$Endo - 1,
  herb_f = demography_climate_distance_flow$Herbivory - 1,
  size_f = demography_climate_distance_flow$log_size_t0,
  y_f = demography_climate_distance_flow$flow_t1,
  n_f = nrow(demography_climate_distance_flow)
)
## Running the stan model
# fit_allsites_flow_aghy_ppt <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/flowering.stan",
#   data = demography_flow_aghy_ppt,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control)

fit_allsites_flow_aghy_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/is4m1b6je38i68dfn6bkf/fit_allsites_flow_aghy_ppt.rds?rlkey=oc3jeswl8jpfkrjb0r20jfnng&dl=1"))

summary(fit_allsites_flow_aghy_ppt)$summary[, c("Rhat", "n_eff")]
posterior_flow_aghy_ppt <- as.array(fit_allsites_flow_aghy_ppt) # Converts to an array
bayesplot::mcmc_trace(posterior_flow_aghy_ppt,
  pars = quote_bare(
    b0_f[1], b0_f[2], b0_f[3],
    bendo_f[1], bendo_f[2], bendo_f[3],
    bherb_f[1], bherb_f[2], bherb_f[3],
    bclim_f[1], bclim_f[2], bclim_f[3],
    bendoclim_f[1], bendoclim_f[2], bendoclim_f[3],
    bendoherb_f[1], bendoherb_f[2], bendoherb_f[3],
    bclim2_f[1], bclim2_f[2], bclim2_f[3],
    bendoclim2_f[1], bendoclim2_f[2], bendoclim2_f[3]
  )
) + theme_bw()

## Compute predictions using posterior draws
posterior_samples_flow_ppt <- rstan::extract(fit_allsites_flow_aghy_ppt)
n_draws <- 1000 # Number of posterior samples to use
pred_data_flow_ppt <- as.data.frame(demography_flow_aghy_ppt)
pred_matrix_flow_ppt <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_flow_ppt))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_flow_ppt)) {
    species <- pred_data_flow_ppt$species_f[j]
    pred_matrix_flow_ppt[i, j] <- posterior_samples_flow_ppt$b0_f[i, species] +
      posterior_samples_flow_ppt$bendo_f[i, species] * pred_data_flow_ppt$endo_f[j] +
      posterior_samples_flow_ppt$bclim_f[i, species] * pred_data_flow_ppt$clim_f[j] +
      posterior_samples_flow_ppt$bherb_f[i, species] * pred_data_flow_ppt$herb_f[j] +
      posterior_samples_flow_ppt$bendoclim_f[i, species] * pred_data_flow_ppt$clim_f[j] * pred_data_flow_ppt$endo_f[j] +
      posterior_samples_flow_ppt$bendoherb_f[i, species] * pred_data_flow_ppt$endo_f[j] * pred_data_flow_ppt$herb_f[j] +
      posterior_samples_flow_ppt$bclim2_f[i, species] * pred_data_flow_ppt$clim_f[j]^2 +
      posterior_samples_flow_ppt$bendoclim2_f[i, species] * pred_data_flow_ppt$endo_f[j] * pred_data_flow_ppt$clim_f[j]^2
  }
}

# Convert to probability scale
pred_prob_flow_ppt <- exp(pred_matrix_flow_ppt)
# Compute mean and 95% credible interval
pred_data_flow_ppt$mean_flowth <- apply(pred_prob_flow_ppt, 2, mean)
pred_data_flow_ppt$lower_ci <- apply(pred_prob_flow_ppt, 2, quantile, probs = 0.025)
pred_data_flow_ppt$upper_ci <- apply(pred_prob_flow_ppt, 2, quantile, probs = 0.975)
# Plot predicted survival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Prflow_ppt.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_flow_ppt, aes(x = exp(clim_f), y = mean_flowth, color = factor(endo_f), fill = factor(endo_f))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_f)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_f), y = y_f, color = factor(endo_f)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  facet_grid(species_f ~ herb_f, labeller = labeller(
    species_f = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_f = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Preciptation (mm)",
    y = " Inflorescence",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.55),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()
# fit_allsites_flow_aghy_pet <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/flowering.stan",
#   data = data_sites_flow_aghy_pet,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control)

fit_allsites_flow_aghy_pet <- readRDS(url("https://www.dropbox.com/scl/fi/8u4tr13xs7jennhsrnbpq/fit_allsites_flow_aghy_pet.rds?rlkey=f2e65i0lu0j81abyiwwtxl2ss&dl=1"))

summary(fit_allsites_flow_aghy_pet)$summary[, c("Rhat", "n_eff")]
posterior_flow_aghy_pet <- as.array(fit_allsites_flow_aghy_pet) # Converts to an array
bayesplot::mcmc_trace(posterior_flow_aghy_pet,
  pars = quote_bare(
    b0_f[1], b0_f[2], b0_f[3],
    bendo_f[1], bendo_f[2], bendo_f[3],
    bherb_f[1], bherb_f[2], bherb_f[3],
    bclim_f[1], bclim_f[2], bclim_f[3],
    bendoclim_f[1], bendoclim_f[2], bendoclim_f[3],
    bendoherb_f[1], bendoherb_f[2], bendoherb_f[3],
    bclim2_f[1], bclim2_f[2], bclim2_f[3],
    bendoclim2_f[1], bendoclim2_f[2], bendoclim2_f[3]
  )
) + theme_bw()

## Compute predictions using posterior draws
posterior_samples_flow_pet <- rstan::extract(fit_allsites_flow_aghy_pet)
n_draws <- 1000 # Number of posterior samples to use
pred_data_flow_pet <- as.data.frame(data_sites_flow_aghy_pet)
pred_matrix_flow_pet <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_flow_pet))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_flow_pet)) {
    species <- pred_data_flow_pet$species_f[j]
    pred_matrix_flow_pet[i, j] <- posterior_samples_flow_pet$b0_f[i, species] +
      posterior_samples_flow_pet$bendo_f[i, species] * pred_data_flow_pet$endo_f[j] +
      posterior_samples_flow_pet$bclim_f[i, species] * pred_data_flow_pet$clim_f[j] +
      posterior_samples_flow_pet$bherb_f[i, species] * pred_data_flow_pet$herb_f[j] +
      posterior_samples_flow_pet$bendoclim_f[i, species] * pred_data_flow_pet$clim_f[j] * pred_data_flow_pet$endo_f[j] +
      posterior_samples_flow_pet$bendoherb_f[i, species] * pred_data_flow_pet$endo_f[j] * pred_data_flow_pet$herb_f[j] +
      posterior_samples_flow_pet$bclim2_f[i, species] * pred_data_flow_pet$clim_f[j]^2 +
      posterior_samples_flow_pet$bendoclim2_f[i, species] * pred_data_flow_pet$endo_f[j] * pred_data_flow_pet$clim_f[j]^2
  }
}

# Convert logits to probability scale
pred_prob_flow_pet <- exp(pred_matrix_flow_pet)
# Compute mean and 95% credible interval
pred_data_flow_pet$mean_flowering <- apply(pred_prob_flow_pet, 2, mean)
pred_data_flow_pet$lower_ci <- apply(pred_prob_flow_pet, 2, quantile, probs = 0.025)
pred_data_flow_pet$upper_ci <- apply(pred_prob_flow_pet, 2, quantile, probs = 0.975)
# Plot predicted flowth  with credible intervals

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Prflow_pet.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_flow_pet, aes(x = exp(clim_f), y = mean_flowering, color = factor(endo_f), fill = factor(endo_f))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_f)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_f), y = y_f, color = factor(endo_f)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  facet_grid(species_f ~ herb_f, labeller = labeller(
    species_f = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_f = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Preciptation (mm)",
    y = "Inflorescence ",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.2, 0.55),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

# fit_allsites_flow_aghy_spei <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/flowering.stan",
#   data = data_sites_flow_aghy_spei,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control)

fit_allsites_flow_aghy_spei <- readRDS(url("https://www.dropbox.com/scl/fi/yaybwm6qn65vg0swkhing/fit_allsites_flow_aghy_spei.rds?rlkey=ncu9cicj1j658r7usohgvxhcw&dl=1"))

summary(fit_allsites_flow_aghy_spei)$summary[, c("Rhat", "n_eff")]
posterior_flow_aghy_spei <- as.array(fit_allsites_flow_aghy_spei) # Converts to an array
bayesplot::mcmc_trace(posterior_flow_aghy_spei,
  pars = quote_bare(
    b0_f[1], b0_f[2], b0_f[3],
    bendo_f[1], bendo_f[2], bendo_f[3],
    bherb_f[1], bherb_f[2], bherb_f[3],
    bclim_f[1], bclim_f[2], bclim_f[3],
    bendoclim_f[1], bendoclim_f[2], bendoclim_f[3],
    bendoherb_f[1], bendoherb_f[2], bendoherb_f[3],
    bclim2_f[1], bclim2_f[2], bclim2_f[3],
    bendoclim2_f[1], bendoclim2_f[2], bendoclim2_f[3]
  )
) + theme_bw()

## Compute predictions using posterior draws
posterior_samples_flow_spei <- rstan::extract(fit_allsites_flow_aghy_spei)
n_draws <- 1000 # Number of posterior samples to use
pred_data_flow_spei <- as.data.frame(data_sites_flow_aghy_spei)
pred_matrix_flow_spei <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_flow_spei))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_flow_spei)) {
    species <- pred_data_flow_spei$species_f[j]
    pred_matrix_flow_spei[i, j] <- posterior_samples_flow_spei$b0_f[i, species] +
      posterior_samples_flow_spei$bendo_f[i, species] * pred_data_flow_spei$endo_f[j] +
      posterior_samples_flow_spei$bclim_f[i, species] * pred_data_flow_spei$clim_f[j] +
      posterior_samples_flow_spei$bherb_f[i, species] * pred_data_flow_spei$herb_f[j] +
      posterior_samples_flow_spei$bendoclim_f[i, species] * pred_data_flow_spei$clim_f[j] * pred_data_flow_spei$endo_f[j] +
      posterior_samples_flow_spei$bendoherb_f[i, species] * pred_data_flow_spei$endo_f[j] * pred_data_flow_spei$herb_f[j] +
      posterior_samples_flow_spei$bclim2_f[i, species] * pred_data_flow_spei$clim_f[j]^2 +
      posterior_samples_flow_spei$bendoclim2_f[i, species] * pred_data_flow_spei$endo_f[j] * pred_data_flow_spei$clim_f[j]^2
  }
}

# Convert logits to probability scale
pred_prob_flow_spei <- exp(pred_matrix_flow_spei)
# Compute mean and 95% credible interval
pred_data_flow_spei$mean_flowering <- apply(pred_prob_flow_spei, 2, mean)
pred_data_flow_spei$lower_ci <- apply(pred_prob_flow_spei, 2, quantile, probs = 0.025)
pred_data_flow_spei$upper_ci <- apply(pred_prob_flow_spei, 2, quantile, probs = 0.975)
# Plot predicted flowvival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Prflow_spei.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_flow_spei, aes(x = exp(clim_f), y = mean_flowering, color = factor(endo_f), fill = factor(endo_f))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_f)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_f), y = y_f, color = factor(endo_f)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  facet_grid(species_f ~ herb_f, labeller = labeller(
    species_f = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_f = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Standardised precipitation-evapotranspiration index",
    y = "Infloresnce",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.9),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

# fit_allsites_flow_aghy_distance <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/flowering.stan",
#   data = data_sites_flow_aghy_distance,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control)

fit_allsites_flow_aghy_distance <- readRDS(url("https://www.dropbox.com/scl/fi/lufyihvzkf7gbl8tfyl1n/fit_allsites_flow_aghy_distance.rds?rlkey=aknvr6uyc5ntyf9m033ycrzrx&dl=1"))

summary(fit_allsites_flow_aghy_distance)$summary[, c("Rhat", "n_eff")]
posterior_flow_aghy_distance <- as.array(fit_allsites_flow_aghy_distance) # Converts to an array
bayesplot::mcmc_trace(posterior_flow_aghy_distance,
  pars = quote_bare(
    b0_f[1], b0_f[2], b0_f[3],
    bendo_f[1], bendo_f[2], bendo_f[3],
    bherb_f[1], bherb_f[2], bherb_f[3],
    bclim_f[1], bclim_f[2], bclim_f[3],
    bendoclim_f[1], bendoclim_f[2], bendoclim_f[3],
    bendoherb_f[1], bendoherb_f[2], bendoherb_f[3],
    bclim2_f[1], bclim2_f[2], bclim2_f[3],
    bendoclim2_f[1], bendoclim2_f[2], bendoclim2_f[3]
  )
) + theme_bw()
## Compute predictions using posterior draws
posterior_samples_flow_distance <- rstan::extract(fit_allsites_flow_aghy_distance)
n_draws <- 1000 # Number of posterior samples to use
pred_data_flow_distance <- as.data.frame(data_sites_flow_aghy_distance)
pred_matrix_flow_distance <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_flow_distance))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_flow_distance)) {
    species <- pred_data_flow_distance$species_f[j]
    pred_matrix_flow_distance[i, j] <- posterior_samples_flow_distance$b0_f[i, species] +
      posterior_samples_flow_distance$bendo_f[i, species] * pred_data_flow_distance$endo_f[j] +
      posterior_samples_flow_distance$bclim_f[i, species] * pred_data_flow_distance$clim_f[j] +
      posterior_samples_flow_distance$bherb_f[i, species] * pred_data_flow_distance$herb_f[j] +
      posterior_samples_flow_distance$bendoclim_f[i, species] * pred_data_flow_distance$clim_f[j] * pred_data_flow_distance$endo_f[j] +
      posterior_samples_flow_distance$bendoherb_f[i, species] * pred_data_flow_distance$endo_f[j] * pred_data_flow_distance$herb_f[j] +
      posterior_samples_flow_distance$bclim2_f[i, species] * pred_data_flow_distance$clim_f[j]^2 +
      posterior_samples_flow_distance$bendoclim2_f[i, species] * pred_data_flow_distance$endo_f[j] * pred_data_flow_distance$clim_f[j]^2
  }
}

# Convert logits to probability scale
pred_prob_flow_distance <- exp(pred_matrix_flow_distance)
# Compute mean and 95% credible interval
pred_data_flow_distance$mean_flowering <- apply(pred_prob_flow_distance, 2, mean)
pred_data_flow_distance$lower_ci <- apply(pred_prob_flow_distance, 2, quantile, probs = 0.025)
pred_data_flow_distance$upper_ci <- apply(pred_prob_flow_distance, 2, quantile, probs = 0.975)
# Plot predicted flowvival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Prflow_distance.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_flow_distance, aes(x = exp(clim_f), y = mean_flowering, color = factor(endo_f), fill = factor(endo_f))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_f)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_f), y = y_f, color = factor(endo_f)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  facet_grid(species_f ~ herb_f, labeller = labeller(
    species_f = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_f = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Mahalanobis distance",
    y = "Infloresnce",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.35, 0.55),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

# fit_allsites_flow_aghy_distance_linear <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/flowering_distance.stan",
#   data = data_sites_flow_aghy_distance,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control)

fit_allsites_flow_aghy_distance_linear <- readRDS(url("https://www.dropbox.com/scl/fi/bxcpb1vacfdr2tq368arv/fit_allsites_flow_aghy_distance_linear.rds?rlkey=13v9pfkt59ybfjlyntm8r54ef&dl=1"))

summary(fit_allsites_flow_aghy_distance_linear)$summary[, c("Rhat", "n_eff")]
posterior_flow_aghy_distance <- as.array(fit_allsites_flow_aghy_distance_linear) # Converts to an array
bayesplot::mcmc_trace(posterior_flow_aghy_distance,
  pars = quote_bare(
    b0_f[1], b0_f[2], b0_f[3],
    bendo_f[1], bendo_f[2], bendo_f[3],
    bherb_f[1], bherb_f[2], bherb_f[3],
    bclim_f[1], bclim_f[2], bclim_f[3],
    bendoclim_f[1], bendoclim_f[2], bendoclim_f[3],
    bendoherb_f[1], bendoherb_f[2], bendoherb_f[3]
  )
) + theme_bw()
## Compute predictions using posterior draws
posterior_samples_flow_distance <- rstan::extract(fit_allsites_flow_aghy_distance_linear)
n_draws <- 1000 # Number of posterior samples to use
pred_data_flow_distance <- as.data.frame(data_sites_flow_aghy_distance)
pred_matrix_flow_distance <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_flow_distance))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_flow_distance)) {
    species <- pred_data_flow_distance$species_f[j]
    pred_matrix_flow_distance[i, j] <- posterior_samples_flow_distance$b0_f[i, species] +
      posterior_samples_flow_distance$bendo_f[i, species] * pred_data_flow_distance$endo_f[j] +
      posterior_samples_flow_distance$bclim_f[i, species] * pred_data_flow_distance$clim_f[j] +
      posterior_samples_flow_distance$bherb_f[i, species] * pred_data_flow_distance$herb_f[j] +
      posterior_samples_flow_distance$bendoclim_f[i, species] * pred_data_flow_distance$clim_f[j] * pred_data_flow_distance$endo_f[j] +
      posterior_samples_flow_distance$bendoherb_f[i, species] * pred_data_flow_distance$endo_f[j] * pred_data_flow_distance$herb_f[j]
  }
}

# Convert logits to probability scale
pred_prob_flow_distance <- pred_matrix_flow_distance
# Compute mean and 95% credible interval
pred_data_flow_distance$mean_flowering <- apply(pred_prob_flow_distance, 2, mean)
pred_data_flow_distance$lower_ci <- apply(pred_prob_flow_distance, 2, quantile, probs = 0.025)
pred_data_flow_distance$upper_ci <- apply(pred_prob_flow_distance, 2, quantile, probs = 0.975)
# Plot predicted flowvival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Prflow_distance_l.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_flow_distance, aes(x = exp(clim_f), y = mean_flowering, color = factor(endo_f), fill = factor(endo_f))) +
  # geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_f)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_f), y = y_f, color = factor(endo_f)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.3) + # Linear regression with confidence interval
  facet_grid(species_f ~ herb_f, labeller = labeller(
    species_f = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_f = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Mahalanobis distance",
    y = "Infloresnce",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.35, 0.55),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()


## Save RDS file for further use
# saveRDS(fit_allsites_flow_aghy_ppt, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_flow_aghy_ppt.rds')
# saveRDS(fit_allsites_flow_aghy_pet, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_flow_aghy_pet.rds')
# saveRDS(fit_allsites_flow_aghy_spei, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_flow_aghy_spei.rds')
# saveRDS(fit_allsites_flow_aghy_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_flow_aghy_distance.rds')
# saveRDS(fit_allsites_flow_aghy_distance_linear, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_flow_aghy_distance_linear.rds')

# Spikelet----
demography_climate_distance %>%
  filter(Species %in% c("ELVI", "POAU")) %>%
  subset(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot, site_plot, Endo, Herbivory,
    tiller_t, spikelet_t1, mean_ppt, mean_pet, mean_spei, distance
  ) %>%
  na.omit() %>%
  mutate(
    Site = Site %>% as.factor() %>% as.numeric(),
    Species = Species %>% as.factor() %>% as.numeric(),
    Plot = Plot %>% as.factor() %>% as.numeric(),
    site_plot = site_plot %>% as.factor() %>% as.numeric(),
    Endo = Endo %>% as.factor() %>% as.numeric(),
    Herbivory = Herbivory %>% as.factor() %>% as.numeric(),
    Population = Population %>% as.factor() %>% as.numeric()
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    spi_t1 = spikelet_t1,
    ppt = log(mean_ppt),
    pet = log(mean_pet),
    spei = mean_spei,
    distance = log(distance)
  ) -> demography_climate_distance_spik

# hist(demography_climate_distance_spike$spikelet_t1,main="")
demography_spik_aghy_ppt <- list(
  n_species = demography_climate_distance_spik$Species %>% n_distinct(),
  n_sites = demography_climate_distance_spik$Site %>% n_distinct(),
  n_pops = demography_climate_distance_spik$Population %>% n_distinct(),
  # flowering data
  n_plot_spk = demography_climate_distance_spik$Plot %>% n_distinct(),
  species_spk = demography_climate_distance_spik$Species,
  site_spk = demography_climate_distance_spik$Site,
  pop_spk = demography_climate_distance_spik$Population,
  plot_spk = demography_climate_distance_spik$Plot,
  clim_spk = as.vector(demography_climate_distance_spik$ppt),
  endo_spk = demography_climate_distance_spik$Endo - 1,
  herb_spk = demography_climate_distance_spik$Herbivory - 1,
  size_spk = demography_climate_distance_spik$log_size_t0,
  y_spk = demography_climate_distance_spik$spikelet_t1,
  n_spk = nrow(demography_climate_distance_spik)
)

demography_spik_aghy_pet <- list(
  n_species = demography_climate_distance_spik$Species %>% n_distinct(),
  n_sites = demography_climate_distance_spik$Site %>% n_distinct(),
  n_pops = demography_climate_distance_spik$Population %>% n_distinct(),
  # flowering data
  n_plot_spk = demography_climate_distance_spik$Plot %>% n_distinct(),
  species_spk = demography_climate_distance_spik$Species,
  site_spk = demography_climate_distance_spik$Site,
  pop_spk = demography_climate_distance_spik$Population,
  plot_spk = demography_climate_distance_spik$Plot,
  clim_spk = as.vector(demography_climate_distance_spik$pet),
  endo_spk = demography_climate_distance_spik$Endo - 1,
  herb_spk = demography_climate_distance_spik$Herbivory - 1,
  size_spk = demography_climate_distance_spik$log_size_t0,
  y_spk = demography_climate_distance_spik$spikelet_t1,
  n_spk = nrow(demography_climate_distance_spik)
)
demography_spik_aghy_spei <- list(
  n_species = demography_climate_distance_spik$Species %>% n_distinct(),
  n_sites = demography_climate_distance_spik$Site %>% n_distinct(),
  n_pops = demography_climate_distance_spik$Population %>% n_distinct(),
  # flowering data
  n_plot_spk = demography_climate_distance_spik$Plot %>% n_distinct(),
  species_spk = demography_climate_distance_spik$Species,
  site_spk = demography_climate_distance_spik$Site,
  pop_spk = demography_climate_distance_spik$Population,
  plot_spk = demography_climate_distance_spik$Plot,
  clim_spk = as.vector(demography_climate_distance_spik$spei),
  endo_spk = demography_climate_distance_spik$Endo - 1,
  herb_spk = demography_climate_distance_spik$Herbivory - 1,
  size_spk = demography_climate_distance_spik$log_size_t0,
  y_spk = demography_climate_distance_spik$spikelet_t1,
  n_spk = nrow(demography_climate_distance_spik)
)
demography_spik_aghy_distance <- list(
  n_species = demography_climate_distance_spik$Species %>% n_distinct(),
  n_sites = demography_climate_distance_spik$Site %>% n_distinct(),
  n_pops = demography_climate_distance_spik$Population %>% n_distinct(),
  # flowering data
  n_plot_spk = demography_climate_distance_spik$Plot %>% n_distinct(),
  species_spk = demography_climate_distance_spik$Species,
  site_spk = demography_climate_distance_spik$Site,
  pop_spk = demography_climate_distance_spik$Population,
  plot_spk = demography_climate_distance_spik$Plot,
  clim_spk = as.vector(demography_climate_distance_spik$distance),
  endo_spk = demography_climate_distance_spik$Endo - 1,
  herb_spk = demography_climate_distance_spik$Herbivory - 1,
  size_spk = demography_climate_distance_spik$log_size_t0,
  y_spk = demography_climate_distance_spik$spikelet_t1,
  n_spk = nrow(demography_climate_distance_spik)
)

mean_value <- mean(demography_spik_aghy_ppt$y_spk, na.rm = TRUE)
variance_value <- var(demography_spik_aghy_ppt$y_spk, na.rm = TRUE)

# fit_allsites_spik_aghy_ppt <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/spikelet.stan",
#   data = demography_spik_aghy_ppt,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control =sim_pars$control)

fit_allsites_spik_aghy_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/8rtir221u5ml997h9usgf/fit_allsites_spik_aghy_ppt.rds?rlkey=4e4ya6tnhnosqqhwu76hfpb2x&dl=1"))

summary(fit_allsites_spik_aghy_ppt)$summary[, c("Rhat", "n_eff")]
posterior_spik_aghy_ppt <- as.array(fit_allsites_spik_aghy_ppt) # Converts to an array
bayesplot::mcmc_trace(posterior_spik_aghy_ppt,
  pars = quote_bare(
    b0_spk[1], b0_spk[2],
    bendo_spk[1], bendo_spk[2],
    bherb_spk[1], bherb_spk[2],
    bclim_spk[1], bclim_spk[2],
    bendoclim_spk[1], bendoclim_spk[2],
    bendoherb_spk[1], bendoherb_spk[2],
    bclim2_spk[1], bclim2_spk[2],
    bendoclim2_spk[1], bendoclim2_spk[2]
  )
) + theme_bw()


## Compute predictions using posterior draws
posterior_samples_spik_ppt <- rstan::extract(fit_allsites_spik_aghy_ppt)
n_draws <- 1000 # Number of posterior samples to use
pred_data_spik_ppt <- as.data.frame(demography_spik_aghy_ppt)
pred_matrix_spik_ppt <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_spik_ppt))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_spik_ppt)) {
    species <- pred_data_spik_ppt$species_spk[j]
    pred_matrix_spik_ppt[i, j] <- posterior_samples_spik_ppt$b0_spk[i, species] +
      posterior_samples_spik_ppt$bendo_spk[i, species] * pred_data_spik_ppt$endo_spk[j] +
      posterior_samples_spik_ppt$bclim_spk[i, species] * pred_data_spik_ppt$clim_spk[j] +
      posterior_samples_spik_ppt$bherb_spk[i, species] * pred_data_spik_ppt$herb_spk[j] +
      posterior_samples_spik_ppt$bendoclim_spk[i, species] * pred_data_spik_ppt$clim_spk[j] * pred_data_spik_ppt$endo_spk[j] +
      posterior_samples_spik_ppt$bendoherb_spk[i, species] * pred_data_spik_ppt$endo_spk[j] * pred_data_spik_ppt$herb_spk[j] +
      posterior_samples_spik_ppt$bclim2_spk[i, species] * pred_data_spik_ppt$clim_spk[j]^2 +
      posterior_samples_spik_ppt$bendoclim2_spk[i, species] * pred_data_spik_ppt$endo_spk[j] * pred_data_spik_ppt$clim_spk[j]^2
  }
}

# Convert to probability scale
pred_prob_spik_ppt <- exp(pred_matrix_spik_ppt)
# Compute mean and 95% credible interval
pred_data_spik_ppt$mean_spikelet <- apply(pred_prob_spik_ppt, 2, mean)
pred_data_spik_ppt$lower_ci <- apply(pred_prob_spik_ppt, 2, quantile, probs = 0.025)
pred_data_spik_ppt$upper_ci <- apply(pred_prob_spik_ppt, 2, quantile, probs = 0.975)
# Plot predicted survival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Prspik_ppt.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_spik_ppt, aes(x = exp(clim_spk), y = mean_spikelet, color = factor(endo_spk), fill = factor(endo_spk))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_spk)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_spk), y = y_spk, color = factor(endo_spk)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  facet_grid(species_spk ~ herb_spk, labeller = labeller(
    species_spk = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_spk = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Preciptation (mm)",
    y = " Spikelets",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.91),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()
# fit_allsites_spik_aghy_pet <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/spikelet.stan",
#   data = demography_spik_aghy_pet,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control =sim_pars$control )

fit_allsites_spik_aghy_pet <- readRDS(url("https://www.dropbox.com/scl/fi/wmn81q56ya2ykf0hk4rrg/fit_allsites_spik_aghy_pet.rds?rlkey=1ifut2cb3zdhh19qm8t53mefh&dl=1"))

summary(fit_allsites_spik_aghy_pet)$summary[, c("Rhat", "n_eff")]
posterior_spik_aghy_pet <- as.array(fit_allsites_spik_aghy_pet) # Converts to an array
bayesplot::mcmc_trace(posterior_spik_aghy_pet,
  pars = quote_bare(
    b0_spk[1], b0_spk[2],
    bendo_spk[1], bendo_spk[2],
    bherb_spk[1], bherb_spk[2],
    bclim_spk[1], bclim_spk[2],
    bendoclim_spk[1], bendoclim_spk[2],
    bendoherb_spk[1], bendoherb_spk[2],
    bclim2_spk[1], bclim2_spk[2],
    bendoclim2_spk[1], bendoclim2_spk[2]
  )
) + theme_bw()

## Compute predictions using posterior draws
posterior_samples_spik_pet <- rstan::extract(fit_allsites_spik_aghy_pet)
n_draws <- 1000 # Number of posterior samples to use
pred_data_spik_pet <- as.data.frame(demography_spik_aghy_pet)
pred_matrix_spik_pet <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_spik_pet))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_spik_pet)) {
    species <- pred_data_spik_pet$species_spk[j]
    pred_matrix_spik_pet[i, j] <- posterior_samples_spik_pet$b0_spk[i, species] +
      posterior_samples_spik_pet$bendo_spk[i, species] * pred_data_spik_pet$endo_spk[j] +
      posterior_samples_spik_pet$bclim_spk[i, species] * pred_data_spik_pet$clim_spk[j] +
      posterior_samples_spik_pet$bherb_spk[i, species] * pred_data_spik_pet$herb_spk[j] +
      posterior_samples_spik_pet$bendoclim_spk[i, species] * pred_data_spik_pet$clim_spk[j] * pred_data_spik_pet$endo_spk[j] +
      posterior_samples_spik_pet$bendoherb_spk[i, species] * pred_data_spik_pet$endo_spk[j] * pred_data_spik_pet$herb_spk[j] +
      posterior_samples_spik_pet$bclim2_spk[i, species] * pred_data_spik_pet$clim_spk[j]^2 +
      posterior_samples_spik_pet$bendoclim2_spk[i, species] * pred_data_spik_pet$endo_spk[j] * pred_data_spik_pet$clim_spk[j]^2
  }
}

# Convert logits to probability scale
pred_prob_spik_pet <- exp(pred_matrix_spik_pet)
# Compute mean and 95% credible interval
pred_data_spik_pet$mean_spikelet <- apply(pred_prob_spik_pet, 2, mean)
pred_data_spik_pet$lower_ci <- apply(pred_prob_spik_pet, 2, quantile, probs = 0.025)
pred_data_spik_pet$upper_ci <- apply(pred_prob_spik_pet, 2, quantile, probs = 0.975)
# Plot predicted spikth  with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Prspik_pet.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_spik_pet, aes(x = exp(clim_spk), y = mean_spikelet, color = factor(endo_spk), fill = factor(endo_spk))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_spk)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_spk), y = y_spk, color = factor(endo_spk)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  facet_grid(species_spk ~ herb_spk, labeller = labeller(
    species_spk = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_spk = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Potential Evapotranspiration",
    y = "Inflorescence ",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.38, 0.40),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()
# fit_allsites_spik_aghy_spei <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/spikelet.stan",
#   data = demography_spik_aghy_spei,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control
# )

fit_allsites_spik_aghy_spei <- readRDS(url("https://www.dropbox.com/scl/fi/swhi510v2mhnlb62xrcvo/fit_allsites_spik_aghy_spei.rds?rlkey=s9szgdfhokb7jjn7zxix5t93e&dl=1"))

summary(fit_allsites_spik_aghy_spei)$summary[, c("Rhat", "n_eff")]
posterior_spik_aghy_spei <- as.array(fit_allsites_spik_aghy_spei) # Converts to an array
bayesplot::mcmc_trace(posterior_spik_aghy_spei,
  pars = quote_bare(
    b0_spk[1], b0_spk[2],
    bendo_spk[1], bendo_spk[2],
    bherb_spk[1], bherb_spk[2],
    bclim_spk[1], bclim_spk[2],
    bendoclim_spk[1], bendoclim_spk[2],
    bendoherb_spk[1], bendoherb_spk[2],
    bclim2_spk[1], bclim2_spk[2],
    bendoclim2_spk[1], bendoclim2_spk[2]
  )
) + theme_bw()

## Compute predictions using posterior draws
posterior_samples_spik_spei <- rstan::extract(fit_allsites_spik_aghy_spei)
n_draws <- 1000 # Number of posterior samples to use
pred_data_spik_spei <- as.data.frame(demography_spik_aghy_spei)
pred_matrix_spik_spei <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_spik_spei))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_spik_spei)) {
    species <- pred_data_spik_spei$species_spk[j]
    pred_matrix_spik_spei[i, j] <- posterior_samples_spik_spei$b0_spk[i, species] +
      posterior_samples_spik_spei$bendo_spk[i, species] * pred_data_spik_spei$endo_spk[j] +
      posterior_samples_spik_spei$bclim_spk[i, species] * pred_data_spik_spei$clim_spk[j] +
      posterior_samples_spik_spei$bherb_spk[i, species] * pred_data_spik_spei$herb_spk[j] +
      posterior_samples_spik_spei$bendoclim_spk[i, species] * pred_data_spik_spei$clim_spk[j] * pred_data_spik_spei$endo_spk[j] +
      posterior_samples_spik_spei$bendoherb_spk[i, species] * pred_data_spik_spei$endo_spk[j] * pred_data_spik_spei$herb_spk[j] +
      posterior_samples_spik_spei$bclim2_spk[i, species] * pred_data_spik_spei$clim_spk[j]^2 +
      posterior_samples_spik_spei$bendoclim2_spk[i, species] * pred_data_spik_spei$endo_spk[j] * pred_data_spik_spei$clim_spk[j]^2
  }
}

# Convert logits to probability scale
pred_prob_spik_spei <- exp(pred_matrix_spik_spei)
# Compute mean and 95% credible interval
pred_data_spik_spei$mean_spikelet <- apply(pred_prob_spik_spei, 2, mean)
pred_data_spik_spei$lower_ci <- apply(pred_prob_spik_spei, 2, quantile, probs = 0.025)
pred_data_spik_spei$upper_ci <- apply(pred_prob_spik_spei, 2, quantile, probs = 0.975)
# Plot predicted spikvival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Prspik_spei.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_spik_spei, aes(x = exp(clim_spk), y = mean_spikelet, color = factor(endo_spk), fill = factor(endo_spk))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_spk)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_spk), y = y_spk, color = factor(endo_spk)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  facet_grid(species_spk ~ herb_spk, labeller = labeller(
    species_spk = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_spk = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Standardised precipitation-evapotranspiration index",
    y = "Predicted survival probability",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.9),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()
# fit_allsites_spik_aghy_distance <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/spikelet.stan",
#   data = demography_spik_aghy_distance,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control
# )

fit_allsites_spik_aghy_distance <- readRDS(url("https://www.dropbox.com/scl/fi/7apjqo8nris1vlgih4enz/fit_allsites_spik_aghy_distance.rds?rlkey=lhcl48ud0uuve6wetkwx4ov14&dl=1"))

summary(fit_allsites_spik_aghy_distance)$summary[, c("Rhat", "n_eff")]
posterior_spik_aghy_distance <- as.array(fit_allsites_spik_aghy_distance) # Converts to an array
bayesplot::mcmc_trace(posterior_spik_aghy_distance,
  pars = quote_bare(
    b0_spk[1], b0_spk[2],
    bendo_spk[1], bendo_spk[2],
    bherb_spk[1], bherb_spk[2],
    bclim_spk[1], bclim_spk[2],
    bendoclim_spk[1], bendoclim_spk[2],
    bendoherb_spk[1], bendoherb_spk[2],
    bclim2_spk[1], bclim2_spk[2],
    bendoclim2_spk[1], bendoclim2_spk[2]
  )
) + theme_bw()
## Compute predictions using posterior draws
posterior_samples_spik_distance <- rstan::extract(fit_allsites_spik_aghy_distance)
n_draws <- 1000 # Number of posterior samples to use
pred_data_spik_distance <- as.data.frame(demography_spik_aghy_distance)
pred_matrix_spik_distance <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_spik_distance))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_spik_distance)) {
    species <- pred_data_spik_distance$species_spk[j]
    pred_matrix_spik_distance[i, j] <- posterior_samples_spik_distance$b0_spk[i, species] +
      posterior_samples_spik_distance$bendo_spk[i, species] * pred_data_spik_distance$endo_spk[j] +
      posterior_samples_spik_distance$bclim_spk[i, species] * pred_data_spik_distance$clim_spk[j] +
      posterior_samples_spik_distance$bherb_spk[i, species] * pred_data_spik_distance$herb_spk[j] +
      posterior_samples_spik_distance$bendoclim_spk[i, species] * pred_data_spik_distance$clim_spk[j] * pred_data_spik_distance$endo_spk[j] +
      posterior_samples_spik_distance$bendoherb_spk[i, species] * pred_data_spik_distance$endo_spk[j] * pred_data_spik_distance$herb_spk[j] +
      posterior_samples_spik_distance$bclim2_spk[i, species] * pred_data_spik_distance$clim_spk[j]^2 +
      posterior_samples_spik_distance$bendoclim2_spk[i, species] * pred_data_spik_distance$endo_spk[j] * pred_data_spik_distance$clim_spk[j]^2
  }
}

# Convert logits to probability scale
pred_prob_spik_distance <- exp(pred_matrix_spik_distance)
# Compute mean and 95% credible interval
pred_data_spik_distance$mean_spikelet <- apply(pred_prob_spik_distance, 2, mean)
pred_data_spik_distance$lower_ci <- apply(pred_prob_spik_distance, 2, quantile, probs = 0.025)
pred_data_spik_distance$upper_ci <- apply(pred_prob_spik_distance, 2, quantile, probs = 0.975)
# Plot predicted spikvival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Prspik_mh.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_spik_distance, aes(x = exp(clim_spk), y = mean_spikelet, color = factor(endo_spk), fill = factor(endo_spk))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_spk)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_spk), y = y_spk, color = factor(endo_spk)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  facet_grid(species_spk ~ herb_spk, labeller = labeller(
    species_spk = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_spk = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Mahalanobis distance",
    y = "Predicted relative spikth",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.35, 0.40),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  ) # Increase facet label size
dev.off()

# fit_allsites_spik_aghy_distance_linear <- stan(
#   file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/spikelet_distance.stan",
#   data = demography_spik_aghy_distance,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control
# )

fit_allsites_spik_aghy_distance_linear <- readRDS(url("https://www.dropbox.com/scl/fi/7apjqo8nris1vlgih4enz/fit_allsites_spik_aghy_distance.rds?rlkey=lhcl48ud0uuve6wetkwx4ov14&dl=1"))

summary(fit_allsites_spik_aghy_distance_linear)$summary[, c("Rhat", "n_eff")]
posterior_spik_aghy_distance <- as.array(fit_allsites_spik_aghy_distance_linear) # Converts to an array
bayesplot::mcmc_trace(posterior_spik_aghy_distance,
  pars = quote_bare(
    b0_spk[1], b0_spk[2],
    bendo_spk[1], bendo_spk[2],
    bherb_spk[1], bherb_spk[2],
    bclim_spk[1], bclim_spk[2],
    bendoclim_spk[1], bendoclim_spk[2],
    bendoherb_spk[1], bendoherb_spk[2]
  )
) + theme_bw()
## Compute predictions using posterior draws
posterior_samples_spik_distance <- rstan::extract(fit_allsites_spik_aghy_distance)
n_draws <- 1000 # Number of posterior samples to use
pred_data_spik_distance <- as.data.frame(demography_spik_aghy_distance)
pred_matrix_spik_distance <- matrix(NA, nrow = n_draws, ncol = nrow(pred_data_spik_distance))

for (i in 1:n_draws) {
  for (j in 1:nrow(pred_data_spik_distance)) {
    species <- pred_data_spik_distance$species_spk[j]
    pred_matrix_spik_distance[i, j] <- posterior_samples_spik_distance$b0_spk[i, species] +
      posterior_samples_spik_distance$bendo_spk[i, species] * pred_data_spik_distance$endo_spk[j] +
      posterior_samples_spik_distance$bclim_spk[i, species] * pred_data_spik_distance$clim_spk[j] +
      posterior_samples_spik_distance$bherb_spk[i, species] * pred_data_spik_distance$herb_spk[j] +
      posterior_samples_spik_distance$bendoclim_spk[i, species] * pred_data_spik_distance$clim_spk[j] * pred_data_spik_distance$endo_spk[j] +
      posterior_samples_spik_distance$bendoherb_spk[i, species] * pred_data_spik_distance$endo_spk[j] * pred_data_spik_distance$herb_spk[j]
  }
}

# Convert logits to probability scale
pred_prob_spik_distance <- pred_matrix_spik_distance
# Compute mean and 95% credible interval
pred_data_spik_distance$mean_spikelet <- apply(pred_prob_spik_distance, 2, mean)
pred_data_spik_distance$lower_ci <- apply(pred_prob_spik_distance, 2, quantile, probs = 0.025)
pred_data_spik_distance$upper_ci <- apply(pred_prob_spik_distance, 2, quantile, probs = 0.975)
# Plot predicted spikvival probabilities with credible intervals
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Prspik_mh_l.pdf", useDingbats = F, height = 7, width = 5)
ggplot(pred_data_spik_distance, aes(x = exp(clim_spk), y = mean_spikelet, color = factor(endo_spk), fill = factor(endo_spk))) +
  # geom_line() + # Mean prediction
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_spk)), alpha = 0.3, color = NA) + # Credible interval
  # geom_point(aes(x = exp(clim_spk), y = y_spk, color = factor(endo_spk)),
  #   position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  # ) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.3) + # Linear regression with confidence interval
  facet_grid(species_spk ~ herb_spk, labeller = labeller(
    species_spk = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_spk = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Mahalanobis distance",
    y = "Predicted relative spikth",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.35, 0.40),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()
## Save RDS file for further use
# saveRDS(fit_allsites_spik_aghy_ppt, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_spik_aghy_ppt.rds')
# saveRDS(fit_allsites_spik_aghy_pet, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_spik_aghy_pet.rds')
# saveRDS(fit_allsites_spik_aghy_spei, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_spik_aghy_spei.rds')
# saveRDS(fit_allsites_spik_aghy_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_spik_aghy_distance.rds')
# saveRDS(fit_allsites_spik_aghy_distance_linear, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_spik_aghy_distance_linear.rds')
