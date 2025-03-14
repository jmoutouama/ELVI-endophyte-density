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
  # scale_color_manual(values = c("red", "blue")) +  # Custom colors for the species
  # scale_fill_manual(values = c("red", "blue")) +
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
  # scale_color_manual(values = c("red", "blue")) +  # Custom colors for the species
  # scale_fill_manual(values = c("red", "blue")) +
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
  # scale_color_manual(values = c("red", "blue")) +  # Custom colors for the species
  # scale_fill_manual(values = c("red", "blue")) +
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
ggplot(pred_data_sur_ppt, aes(x = exp(clim_s), y = mean_survival, color = factor(endo_s), fill = factor(endo_s))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci,fill=factor(endo_s)), alpha = 0.3,color = NA) +  # Credible interval
  # geom_point(aes(x = exp(clim_s), y = y_s, color = factor(endo_s)),
  #            position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5) +
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
  scale_color_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.8),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  ) # Increase facet label size)

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
ggplot(pred_data_sur_pet, aes(x = exp(clim_s), y = mean_survival, color = factor(endo_s), fill = factor(endo_s))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +  # Credible interval
  # geom_point(aes(x = exp(clim_s), y = y_s, color = factor(endo_s)),
  #            position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5) +
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
  scale_color_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.8),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  ) # Increase facet label size)


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
ggplot(pred_data_sur_spei, aes(x = exp(clim_s), y = mean_survival, color = factor(endo_s), fill = factor(endo_s))) +
  geom_line() + # Mean prediction
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +  # Credible interval
  # geom_point(aes(x = exp(clim_s), y = y_s, color = factor(endo_s)),
  #            position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5) +
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
  scale_color_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.8),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  ) # Increase facet label size)

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
  scale_color_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.8),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  ) # Increase facet label size



## Save RDS file for further use
# saveRDS(fit_allsites_surv_aghy_ppt, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_surv_aghy_ppt.rds')
# saveRDS(fit_allsites_surv_aghy_pet, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_surv_aghy_pet.rds')
# saveRDS(fit_allsites_surv_aghy_spei, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_surv_aghy_spei.rds')
# saveRDS(fit_allsites_surv_aghy_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_surv_aghy_distance.rds')

## Growth----
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
n_draws <- 2000 # Number of posterior samples to use
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
ggplot(pred_data_grow_ppt, aes(x = exp(clim_g), y = mean_growth, color = factor(endo_g), fill = factor(endo_g))) +
  geom_line() + # Mean prediction
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_g)), alpha = 0.3, color = NA) + # Credible interval
  geom_point(aes(x = exp(clim_g), y = y_g, color = factor(endo_g)),
    position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  ) +
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
    legend.position = c(0.9, 0.93),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )


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
ggplot(pred_data_grow_pet, aes(x = exp(clim_g), y = mean_growth, color = factor(endo_g), fill = factor(endo_g))) +
  geom_line() + # Mean prediction
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_g)), alpha = 0.3, color = NA) + # Credible interval
  geom_point(aes(x = exp(clim_g), y = y_g, color = factor(endo_g)),
    position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  ) +
  facet_grid(species_g ~ herb_g, labeller = labeller(
    species_g = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_g = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Preciptation (mm)",
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
  ) # Increase facet label size)

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
ggplot(pred_data_grow_spei, aes(x = exp(clim_g), y = mean_growth, color = factor(endo_g), fill = factor(endo_g))) +
  geom_line() + # Mean prediction
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_g)), alpha = 0.3, color = NA) + # Credible interval
  geom_point(aes(x = exp(clim_g), y = y_g, color = factor(endo_g)),
    position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  ) +
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
ggplot(pred_data_grow_distance, aes(x = exp(clim_g), y = mean_growth, color = factor(endo_g), fill = factor(endo_g))) +
  geom_line() + # Mean prediction
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_g)), alpha = 0.3, color = NA) + # Credible interval
  geom_point(aes(x = exp(clim_g), y = y_g, color = factor(endo_g)),
    position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  ) +
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
  ) # Increase facet label size

## Save RDS file for further use
# saveRDS(fit_allsites_grow_aghy_ppt, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_grow_aghy_ppt.rds')
# saveRDS(fit_allsites_grow_aghy_pet, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_grow_aghy_pet.rds')
# saveRDS(fit_allsites_grow_aghy_spei, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_grow_aghy_spei.rds')
# saveRDS(fit_allsites_grow_aghy_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_grow_aghy_distance.rds')

## Flowering----
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
ggplot(pred_data_flow_ppt, aes(x = exp(clim_f), y = mean_flowth, color = factor(endo_f), fill = factor(endo_f))) +
  geom_line() + # Mean prediction
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_f)), alpha = 0.3, color = NA) + # Credible interval
  geom_point(aes(x = exp(clim_f), y = y_f, color = factor(endo_f)),
    position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  ) +
  facet_grid(species_f ~ herb_f, labeller = labeller(
    species_f = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_f = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Preciptation (mm)",
    y = " Inflorescences",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.55),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )

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
ggplot(pred_data_flow_pet, aes(x = exp(clim_f), y = mean_flowering, color = factor(endo_f), fill = factor(endo_f))) +
  geom_line() + # Mean prediction
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_f)), alpha = 0.3, color = NA) + # Credible interval
  geom_point(aes(x = exp(clim_f), y = y_f, color = factor(endo_f)),
    position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  ) +
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
  scale_color_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.2, 0.55),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  ) # Increase facet label size)


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
ggplot(pred_data_flow_spei, aes(x = exp(clim_f), y = mean_flowering, color = factor(endo_f), fill = factor(endo_f))) +
  geom_line() + # Mean prediction
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_f)), alpha = 0.3, color = NA) + # Credible interval
  geom_point(aes(x = exp(clim_f), y = y_f, color = factor(endo_f)),
    position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  ) +
  facet_grid(species_f ~ herb_f, labeller = labeller(
    species_f = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_f = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Standardised precipitation-evapotranspiration index",
    y = "Predicted survival probability",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.9),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )

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
ggplot(pred_data_flow_distance, aes(x = exp(clim_f), y = mean_flowering, color = factor(endo_f), fill = factor(endo_f))) +
  geom_line() + # Mean prediction
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_f)), alpha = 0.3, color = NA) + # Credible interval
  geom_point(aes(x = exp(clim_f), y = y_f, color = factor(endo_f)),
    position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  ) +
  facet_grid(species_f ~ herb_f, labeller = labeller(
    species_f = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_f = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Mahalanobis distance",
    y = "Predicted relative flowth",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.35, 0.55),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  ) # Increase facet label size



## Save RDS file for further use
# saveRDS(fit_allsites_flow_aghy_ppt, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_flow_aghy_ppt.rds')
# saveRDS(fit_allsites_flow_aghy_pet, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_flow_aghy_pet.rds')
# saveRDS(fit_allsites_flow_aghy_spei, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_flow_aghy_spei.rds')
# saveRDS(fit_allsites_flow_aghy_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_flow_aghy_distance.rds')

## Spikelet----
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

pairs(fit_allsites_spik_aghy_ppt, pars = c("b0_spk", "bclim_spk", "phi_spk"))

summary(fit_allsites_spik_aghy_ppt)$summary[, c("Rhat", "n_eff")]
posterior_spik_aghy_ppt <- as.array(fit_allsites_spik_aghy_ppt) # Converts to an array
bayesplot::mcmc_trace(posterior_spik_aghy_ppt,
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
ggplot(pred_data_spik_ppt, aes(x = exp(clim_spk), y = mean_spikelet, color = factor(endo_spk), fill = factor(endo_spk))) +
  geom_line() + # Mean prediction
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_spk)), alpha = 0.3, color = NA) + # Credible interval
  geom_point(aes(x = exp(clim_spk), y = y_spk, color = factor(endo_spk)),
    position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  ) +
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
  scale_color_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.91),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )

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
pred_data_spik_pet <- as.data.frame(data_sites_spik_aghy_pet)
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
ggplot(pred_data_spik_pet, aes(x = exp(clim_spk), y = mean_spikelet, color = factor(endo_spk), fill = factor(endo_spk))) +
  geom_line() + # Mean prediction
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_spk)), alpha = 0.3, color = NA) + # Credible interval
  geom_point(aes(x = exp(clim_spk), y = y_spk, color = factor(endo_spk)),
    position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  ) +
  facet_grid(species_spk ~ herb_spk, labeller = labeller(
    species_spk = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
    herb_spk = c("0" = "Unfenced", "1" = "Fenced")
  )) + # Facet by species & herbivory
  labs(
    x = "Preciptation (mm)",
    y = "Inflorescence ",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.38, 0.40),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  ) # Increase facet label size)


fit_allsites_spik_aghy_spei <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/spikelet.stan",
  data = data_sites_spik_aghy_spei,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control
)

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
pred_data_spik_spei <- as.data.frame(data_sites_spik_aghy_spei)
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
ggplot(pred_data_spik_spei, aes(x = exp(clim_spk), y = mean_spikelet, color = factor(endo_spk), fill = factor(endo_spk))) +
  geom_line() + # Mean prediction
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_spk)), alpha = 0.3, color = NA) + # Credible interval
  geom_point(aes(x = exp(clim_spk), y = y_spk, color = factor(endo_spk)),
    position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  ) +
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
  scale_color_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.9),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )

fit_allsites_spik_aghy_distance <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/stan/spikelet.stan",
  data = data_sites_spik_aghy_distance,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control
)

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
pred_data_spik_distance <- as.data.frame(data_sites_spik_aghy_distance)
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
ggplot(pred_data_spik_distance, aes(x = exp(clim_spk), y = mean_spikelet, color = factor(endo_spk), fill = factor(endo_spk))) +
  geom_line() + # Mean prediction
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(endo_spk)), alpha = 0.3, color = NA) + # Credible interval
  geom_point(aes(x = exp(clim_spk), y = y_spk, color = factor(endo_spk)),
    position = position_jitter(width = 0.1, height = 0.02), alpha = 0.5
  ) +
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
  scale_color_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "blue", "1" = "red"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.35, 0.40),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  ) # Increase facet label size


## Save RDS file for further use
# saveRDS(fit_allsites_spik_aghy_ppt, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_spik_aghy_ppt.rds')
# saveRDS(fit_allsites_spik_aghy_pet, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_spik_aghy_pet.rds')
# saveRDS(fit_allsites_spik_aghy_spei, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_spik_aghy_spei.rds')
# saveRDS(fit_allsites_spik_aghy_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Endo Model output/fit_allsites_spik_aghy_distance.rds')

## Posterior predictive check for  all models
predS_Tmean <- rstan::extract(fit_allsites_surv_temp_mean, pars = c("predS"))$predS
predS_Tcv <- rstan::extract(fit_allsites_surv_temp_cv, pars = c("predS"))$predS
predS_Wmean <- rstan::extract(fit_allsites_surv_water_mean, pars = c("predS"))$predS
predS_Wcv <- rstan::extract(fit_allsites_surv_water_cv, pars = c("predS"))$predS

predG_Tmean <- rstan::extract(fit_allsites_grow_temp_mean, pars = c("predG"))$predG
sigma_Tmean <- rstan::extract(fit_allsites_grow_temp_mean, pars = c("sigma"))$sigma
predG_Tcv <- rstan::extract(fit_allsites_grow_temp_cv, pars = c("predG"))$predG
sigma_Tcv <- rstan::extract(fit_allsites_grow_temp_cv, pars = c("sigma"))$sigma
predG_Wmean <- rstan::extract(fit_allsites_grow_water_mean, pars = c("predG"))$predG
sigma_Wmean <- rstan::extract(fit_allsites_grow_water_mean, pars = c("sigma"))$sigma
predG_Wcv <- rstan::extract(fit_allsites_grow_water_cv, pars = c("predG"))$predG
sigma_Wcv <- rstan::extract(fit_allsites_grow_water_cv, pars = c("sigma"))$sigma

predF_Tmean <- rstan::extract(fit_allsites_flow_temp_mean, pars = c("predF"))$predF
phi_f_Tmean <- rstan::extract(fit_allsites_flow_temp_mean, pars = c("phi_f"))$phi_f
predF_Tcv <- rstan::extract(fit_allsites_flow_temp_cv, pars = c("predF"))$predF
phi_f_Tcv <- rstan::extract(fit_allsites_flow_temp_cv, pars = c("phi_f"))$phi_f
predF_Wmean <- rstan::extract(fit_allsites_flow_water_mean, pars = c("predF"))$predF
phi_f_Wmean <- rstan::extract(fit_allsites_flow_water_mean, pars = c("phi_f"))$phi_f
predF_Wcv <- rstan::extract(fit_allsites_flow_water_cv, pars = c("predF"))$predF
phi_f_Wcv <- rstan::extract(fit_allsites_flow_water_cv, pars = c("phi_f"))$phi_f

predSP_Tmean <- rstan::extract(fit_allsites_spi_temp_mean, pars = c("predF"))$predF
phi_spk_Tmean <- rstan::extract(fit_allsites_spi_temp_mean, pars = c("phi_spk"))$phi_spk
predSP_Tcv <- rstan::extract(fit_allsites_spi_temp_cv, pars = c("predF"))$predF
phi_spk_Tcv <- rstan::extract(fit_allsites_spi_temp_cv, pars = c("phi_spk"))$phi_spk
predSP_Wmean <- rstan::extract(fit_allsites_spi_water_mean, pars = c("predF"))$predF
phi_spk_Wmean <- rstan::extract(fit_allsites_spi_water_mean, pars = c("phi_spk"))$phi_spk
predSP_Wcv <- rstan::extract(fit_allsites_spi_water_cv, pars = c("predF"))$predF
phi_spk_Wmean <- rstan::extract(fit_allsites_spi_water_cv, pars = c("phi_spk"))$phi_spk

# draw 500 random samples from the joint posterior
n_post_draws <- 500
post_draws <- sample.int(dim(predS_Tmean)[1], n_post_draws)
# set up simulation output
y_s_tm_sim <- matrix(NA, n_post_draws, length(data_sites_surv_temp_mean$y_s))
y_s_tcv_sim <- matrix(NA, n_post_draws, length(data_sites_surv_temp_cv$y_s))
y_s_wm_sim <- matrix(NA, n_post_draws, length(data_sites_surv_water_mean$y_s))
y_s_wcv_sim <- matrix(NA, n_post_draws, length(data_sites_surv_water_cv$y_s))

y_g_tm_sim <- matrix(NA, n_post_draws, length(data_sites_grow_temp_mean$y_g))
y_g_tcv_sim <- matrix(NA, n_post_draws, length(data_sites_grow_temp_cv$y_g))
y_g_wm_sim <- matrix(NA, n_post_draws, length(data_sites_grow_water_mean$y_g))
y_g_wcv_sim <- matrix(NA, n_post_draws, length(data_sites_grow_water_cv$y_g))

y_f_tm_sim <- matrix(NA, n_post_draws, length(data_sites_flow_temp_mean$y_f))
y_f_tcv_sim <- matrix(NA, n_post_draws, length(data_sites_flow_temp_cv$y_f))
y_f_wm_sim <- matrix(NA, n_post_draws, length(data_sites_flow_water_mean$y_f))
y_f_wcv_sim <- matrix(NA, n_post_draws, length(data_sites_flow_water_cv$y_f))

y_sp_tm_sim <- matrix(NA, n_post_draws, length(data_sites_spi_temp_mean$y_spk))
y_sp_tcv_sim <- matrix(NA, n_post_draws, length(data_sites_spi_temp_cv$y_spk))
y_sp_wm_sim <- matrix(NA, n_post_draws, length(data_sites_spi_water_mean$y_spk))
y_sp_wcv_sim <- matrix(NA, n_post_draws, length(data_sites_spi_water_cv$y_spk))


# loop over the posterior and generate new observations
for (i in 1:n_post_draws) {
  print(i)
  ## sample survival data (bernoulli)
  y_s_tm_sim[i, ] <- rbinom(n = length(data_sites_surv_temp_mean$y_s), size = 1, prob = invlogit(predS_Tmean[i, ]))
  y_s_tcv_sim[i, ] <- rbinom(n = length(data_sites_surv_temp_cv$y_s), size = 1, prob = invlogit(predS_Tcv[i, ]))
  y_s_wm_sim[i, ] <- rbinom(n = length(data_sites_surv_water_mean$y_s), size = 1, prob = invlogit(predS_Wmean[i, ]))
  y_s_wcv_sim[i, ] <- rbinom(n = length(data_sites_surv_water_cv$y_s), size = 1, prob = invlogit(predS_Wcv[i, ]))
  ## sample growth data (normal)
  y_g_tm_sim[i, ] <- rnorm(n = length(data_sites_grow_temp_mean$y_g), mean = predG_Tmean[i, ], sd = sigma_Tmean[i])
  y_g_tcv_sim[i, ] <- rnorm(n = length(data_sites_grow_temp_cv$y_g), mean = predG_Tcv[i, ], sd = sigma_Tcv[i])
  y_g_wm_sim[i, ] <- rnorm(n = length(data_sites_grow_water_mean$y_g), mean = predG_Wmean[i, ], sd = sigma_Wmean[i])
  y_g_wcv_sim[i, ] <- rnorm(n = length(data_sites_grow_water_cv$y_g), mean = predG_Wcv[i, ], sd = sigma_Wcv[i])
  ## sample flowering data (negative binomial)
  y_f_tm_sim[i, ] <- rnbinom(n = length(data_sites_flow_temp_mean$y_f), mu = exp(predF_Tmean[i, ]), size = phi_f_Tmean[i])
  y_f_tcv_sim[i, ] <- rnbinom(n = length(data_sites_flow_temp_cv$y_f), mu = exp(predF_Tcv[i, ]), size = phi_f_Tcv[i])
  y_f_wm_sim[i, ] <- rnbinom(n = length(data_sites_flow_water_mean$y_f), mu = exp(predF_Wmean[i, ]), size = phi_f_Wmean[i])
  y_f_wcv_sim[i, ] <- rnbinom(n = length(data_sites_flow_water_cv$y_f), mu = exp(predF_Wcv[i, ]), size = phi_f_Tcv[i])
  ## sample spikelet data (negative binomial)
  y_sp_tm_sim[i, ] <- rnbinom(n = length(data_sites_spi_temp_mean$y_spk), mu = exp(predSP_Tmean[i, ]), size = phi_spk_Tmean[i])
  y_sp_tcv_sim[i, ] <- rnbinom(n = length(data_sites_spi_temp_cv$y_spk), mu = exp(predSP_Tcv[i, ]), size = phi_spk_Tcv)
  y_sp_wm_sim[i, ] <- rnbinom(n = length(data_sites_spi_water_mean$y_spk), mu = exp(predSP_Wmean[i, ]), size = phi_spk_Wmean)
  y_sp_wcv_sim[i, ] <- rnbinom(n = length(data_sites_spi_water_cv$y_spk), mu = exp(predSP_Wcv[i, ]), size = phi_spk_Tcv)
}

# plot ppc overlay
# color_scheme_set("red")
color_scheme_set("blue")
bayesplot::ppc_dens_overlay(data_sites_surv_temp_mean$y_s, y_s_tm_sim) +
  xlab("Survival ") +
  ylab("Density") +
  ggtitle(("Mean temperature")) +
  theme_bw() +
  theme(legend.position = c(0.15, 0.7)) -> ppc_surv_temp_mean
ppc_dens_overlay(data_sites_surv_temp_mean$y_s, y_s_tm_sim) +
  xlab("Survival ") +
  ylab("Density") +
  ggtitle(("CV temperature")) +
  theme_bw() +
  theme(legend.position = "none") -> ppc_surv_temp_cv
ppc_dens_overlay(data_sites_surv_temp_mean$y_s, y_s_tm_sim) +
  xlab("Survival ") +
  ylab("Density") +
  ggtitle(("Mean soil moisture")) +
  theme_bw() +
  theme(legend.position = "none") -> ppc_surv_water_mean
ppc_dens_overlay(data_sites_surv_temp_mean$y_s, y_s_tm_sim) +
  xlab("Survival ") +
  ylab("Density") +
  ggtitle(("CV soil moisture")) +
  theme_bw() +
  theme(legend.position = "none") -> ppc_surv_water_cv
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PPC_survival.pdf", useDingbats = F, height = 5, width = 6)
multiplot(ppc_surv_temp_mean, ppc_surv_temp_cv, ppc_surv_water_mean,
  ppc_surv_water_cv,
  cols = 2
)
dev.off()

bayesplot::ppc_dens_overlay(data_sites_grow_temp_mean$y_g, y_g_tm_sim) +
  xlab("Growth ") +
  ylab("Density") +
  ggtitle(("Mean temperature")) +
  theme_bw() +
  theme(legend.position = c(0.15, 0.7)) -> ppc_grow_temp_mean
ppc_dens_overlay(data_sites_grow_temp_mean$y_g, y_g_tm_sim) +
  xlab("Growth ") +
  ylab("Density") +
  ggtitle(("CV temperature")) +
  theme_bw() +
  theme(legend.position = "none") -> ppc_grow_temp_cv
ppc_dens_overlay(data_sites_grow_temp_mean$y_g, y_g_tm_sim) +
  xlab("Growth ") +
  ylab("Density") +
  ggtitle(("Mean soil moisture")) +
  theme_bw() +
  theme(legend.position = "none") -> ppc_grow_water_mean
ppc_dens_overlay(data_sites_grow_temp_mean$y_g, y_g_tm_sim) +
  xlab("Growth ") +
  ylab("Density") +
  ggtitle(("CV soil moisture")) +
  theme_bw() +
  theme(legend.position = "none") -> ppc_grow_water_cv
# pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PPC_growth.pdf",useDingbats = F,height=5,width=6)
# multiplot(ppc_grow_temp_mean,ppc_grow_temp_cv,ppc_grow_water_mean,
#           ppc_grow_water_cv,cols=2)
# dev.off()

bayesplot::ppc_dens_overlay(data_sites_flow_temp_mean$y_f, y_f_tm_sim) +
  xlab("Flowering") +
  ylab("Density") +
  xlim(-10, 40) +
  ggtitle(("Mean temperature")) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.7)) -> ppc_flow_temp_mean
ppc_dens_overlay(data_sites_flow_temp_mean$y_f, y_f_tm_sim) +
  xlab("Flowering") +
  ylab("Density") +
  xlim(-10, 40) +
  ggtitle(("CV temperature")) +
  theme_bw() +
  theme(legend.position = "none") -> ppc_flow_temp_cv
ppc_dens_overlay(data_sites_flow_temp_mean$y_f, y_f_tm_sim) +
  xlab("Flowering") +
  ylab("Density") +
  xlim(-10, 40) +
  ggtitle(("Mean soil moisture")) +
  theme_bw() +
  theme(legend.position = "none") -> ppc_flow_water_mean
ppc_dens_overlay(data_sites_flow_temp_mean$y_f, y_f_tm_sim) +
  xlab("Flowering") +
  ylab("Density") +
  xlim(-10, 40) +
  ggtitle(("CV soil moisture")) +
  theme_bw() +
  theme(legend.position = "none") -> ppc_flow_water_cv
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PPC_flowering.pdf", useDingbats = F, height = 5, width = 6)
multiplot(ppc_flow_temp_mean, ppc_flow_temp_cv, ppc_flow_water_mean,
  ppc_flow_water_cv,
  cols = 2
)
dev.off()

bayesplot::ppc_dens_overlay(data_sites_spi_temp_mean$y_spk, y_sp_tm_sim) +
  xlab("Spikelet") +
  ylab("Density") +
  xlim(-10, 200) +
  ggtitle(("Mean temperature")) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.7)) -> ppc_spi_temp_mean
ppc_dens_overlay(data_sites_spi_temp_mean$y_spk, y_sp_tm_sim) +
  xlab("Spikelet") +
  ylab("Density") +
  xlim(-10, 200) +
  ggtitle(("CV temperature")) +
  theme_bw() +
  theme(legend.position = "none") -> ppc_spi_temp_cv
ppc_dens_overlay(data_sites_spi_temp_mean$y_spk, y_sp_tm_sim) +
  xlab("Spikelet") +
  ylab("Density") +
  xlim(-10, 200) +
  ggtitle(("Mean soil moisture")) +
  theme_bw() +
  theme(legend.position = "none") -> ppc_spi_water_mean
ppc_dens_overlay(data_sites_spi_temp_mean$y_spk, y_sp_tm_sim) +
  xlab("Spikelet") +
  ylab("Density") +
  xlim(-10, 200) +
  ggtitle(("CV soil moisture")) +
  theme_bw() +
  theme(legend.position = "none") -> ppc_spi_water_cv
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PPC_spikelet.pdf", useDingbats = F, height = 5, width = 6)
multiplot(ppc_spi_temp_mean, ppc_spi_temp_cv, ppc_spi_water_mean,
  ppc_spi_water_cv,
  cols = 2
)
dev.off()


## Model selection -----
## Model comparison based on epdl/looic
### Survival
log_lik_surv_temp_mean <- loo::extract_log_lik(fit_allsites_surv_temp_mean, merge_chains = FALSE)
r_eff_surv_temp_mean <- loo::relative_eff(exp(log_lik_surv_temp_mean))
loo_surv_temp_mean <- loo(log_lik_surv_temp_mean, r_eff = r_eff_surv_temp_mean, cores = 4)
# plot(loo_surv_temp_mean)

log_lik_surv_temp_cv <- loo::extract_log_lik(fit_allsites_surv_temp_cv, merge_chains = FALSE)
r_eff_surv_temp_cv <- loo::relative_eff(exp(log_lik_surv_temp_cv))
loo_surv_temp_cv <- loo(log_lik_surv_temp_cv, r_eff = r_eff_surv_temp_cv, cores = 4)
# plot(loo_surv_temp_cv)

log_lik_surv_water_mean <- loo::extract_log_lik(fit_allsites_surv_water_mean, merge_chains = FALSE)
r_eff_surv_water_mean <- loo::relative_eff(exp(log_lik_surv_water_mean))
loo_surv_water_mean <- loo(log_lik_surv_water_mean, r_eff = r_eff_surv_water_mean, cores = 4)
# plot(loo_surv_water_mean)

log_lik_surv_water_cv <- loo::extract_log_lik(fit_allsites_surv_water_cv, merge_chains = FALSE)
r_eff_surv_water_cv <- loo::relative_eff(exp(log_lik_surv_water_cv))
loo_surv_water_cv <- loo(log_lik_surv_water_cv, r_eff = r_eff_surv_water_cv, cores = 4)
# plot(loo_surv_water_cv)

(comp_surv <- loo::loo_compare(loo_surv_temp_mean, loo_surv_temp_cv, loo_surv_water_mean, loo_surv_water_cv))

### Growth
log_lik_grow_temp_mean <- loo::extract_log_lik(fit_allsites_grow_temp_mean, merge_chains = FALSE)
r_eff_grow_temp_mean <- loo::relative_eff(exp(log_lik_grow_temp_mean))
loo_grow_temp_mean <- loo(log_lik_grow_temp_mean, r_eff = r_eff_grow_temp_mean, cores = 4)
# plot(loo_grow_temp_mean)

log_lik_grow_temp_cv <- loo::extract_log_lik(fit_allsites_grow_temp_cv, merge_chains = FALSE)
r_eff_grow_temp_cv <- loo::relative_eff(exp(log_lik_grow_temp_cv))
loo_grow_temp_cv <- loo(log_lik_grow_temp_cv, r_eff = r_eff_grow_temp_cv, cores = 4)
# plot(loo_grow_temp_cv)

log_lik_grow_water_mean <- loo::extract_log_lik(fit_allsites_grow_water_mean, merge_chains = FALSE)
r_eff_grow_water_mean <- loo::relative_eff(exp(log_lik_grow_water_mean))
loo_grow_water_mean <- loo(log_lik_grow_water_mean, r_eff = r_eff_grow_water_mean, cores = 4)
# plot(loo_grow_water_mean)

log_lik_grow_water_cv <- loo::extract_log_lik(fit_allsites_grow_water_cv, merge_chains = FALSE)
r_eff_grow_water_cv <- loo::relative_eff(exp(log_lik_grow_water_cv))
loo_grow_water_cv <- loo(log_lik_grow_water_cv, r_eff = r_eff_grow_water_cv, cores = 4)
# plot(loo_grow_water_cv)

(comp_grow <- loo::loo_compare(loo_grow_temp_mean, loo_grow_temp_cv, loo_grow_water_mean, loo_grow_water_cv))

## Flowering
log_lik_flow_temp_mean <- loo::extract_log_lik(fit_allsites_flow_temp_mean, merge_chains = FALSE)
r_eff_flow_temp_mean <- loo::relative_eff(exp(log_lik_flow_temp_mean))
loo_flow_temp_mean <- loo(log_lik_flow_temp_mean, r_eff = r_eff_flow_temp_mean, cores = 4)
# plot(loo_flow_temp_mean)

log_lik_flow_temp_cv <- loo::extract_log_lik(fit_allsites_flow_temp_cv, merge_chains = FALSE)
r_eff_flow_temp_cv <- loo::relative_eff(exp(log_lik_flow_temp_cv))
loo_flow_temp_cv <- loo(log_lik_flow_temp_cv, r_eff = r_eff_flow_temp_cv, cores = 4)
# plot(loo_flow_temp_cv)

log_lik_flow_water_mean <- loo::extract_log_lik(fit_allsites_flow_water_mean, merge_chains = FALSE)
r_eff_flow_water_mean <- loo::relative_eff(exp(log_lik_flow_water_mean))
loo_flow_water_mean <- loo(log_lik_flow_water_mean, r_eff = r_eff_flow_water_mean, cores = 4)
# plot(loo_flow_water_mean)

log_lik_flow_water_cv <- loo::extract_log_lik(fit_allsites_flow_water_cv, merge_chains = FALSE)
r_eff_flow_water_cv <- loo::relative_eff(exp(log_lik_flow_water_cv))
loo_flow_water_cv <- loo(log_lik_flow_water_cv, r_eff = r_eff_flow_water_cv, cores = 4)
# plot(loo_flow_water_cv)

(comp_flow <- loo::loo_compare(loo_flow_temp_mean, loo_flow_temp_cv, loo_flow_water_mean, loo_flow_water_cv))

## Spiekelet
log_lik_spi_temp_mean <- loo::extract_log_lik(fit_allsites_spi_temp_mean, merge_chains = FALSE)
r_eff_spi_temp_mean <- loo::relative_eff(exp(log_lik_spi_temp_mean))
loo_spi_temp_mean <- loo(log_lik_spi_temp_mean, r_eff = r_eff_spi_temp_mean, cores = 4)
# plot(loo_spi_temp_mean)

log_lik_spi_temp_cv <- loo::extract_log_lik(fit_allsites_spi_temp_cv, merge_chains = FALSE)
r_eff_spi_temp_cv <- loo::relative_eff(exp(log_lik_spi_temp_cv))
loo_spi_temp_cv <- loo(log_lik_spi_temp_cv, r_eff = r_eff_spi_temp_cv, cores = 4)
# plot(loo_spi_temp_cv)

log_lik_spi_water_mean <- loo::extract_log_lik(fit_allsites_spi_water_mean, merge_chains = FALSE)
r_eff_spi_water_mean <- loo::relative_eff(exp(log_lik_spi_water_mean))
loo_spi_water_mean <- loo(log_lik_spi_water_mean, r_eff = r_eff_spi_water_mean, cores = 4)
# plot(loo_spi_water_mean)

log_lik_spi_water_cv <- loo::extract_log_lik(fit_allsites_spi_water_cv, merge_chains = FALSE)
r_eff_spi_water_cv <- loo::relative_eff(exp(log_lik_spi_water_cv))
loo_spi_water_cv <- loo(log_lik_spi_water_cv, r_eff = r_eff_spi_water_cv, cores = 4)
# plot(loo_spi_water_cv)

(comp_spi <- loo::loo_compare(loo_spi_temp_mean, loo_spi_temp_cv, loo_spi_water_mean, loo_spi_water_cv))

## # Posterior mean values for each vital rate----
## Temp CV
posterior_surv <- as.array(fit_allsites_surv_temp_cv)
# color_scheme_set("red")
surv <- mcmc_intervals(posterior_surv, pars = quote_bare(
  b0_s, bendo_s, bherb_s, btemp_s,
  bendotemp_s, bherbtemp_s, bendoherb_s
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_s", "bendo_s", "bherb_s", "btemp_s",
      "bendotemp_s", "bherbtemp_s", "bendoherb_s"
    ),
    labels = c(
      "b0_s" = "Grand Mean",
      "bendo_s" = "Endophyte",
      "bherb_s" = "Herbivory",
      "bsizesex_s" = "size:sex",
      "btemp_s" = "Temperature",
      "bendotemp_s" = "Endo:Temp",
      "bherbtemp_s" = "Herb:Temp",
      "bendoherb_s" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Survival)") +
  xlim(-8, 8) +
  ggtitle("A") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

posterior_grow <- as.array(fit_allsites_grow_temp_cv)
grow <- mcmc_intervals(posterior_grow, pars = quote_bare(
  b0_g, bendo_g, bherb_g, btemp_g,
  bendotemp_g, bherbtemp_g, bendoherb_g
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_g", "bendo_g", "bherb_g", "btemp_g",
      "bendotemp_g", "bherbtemp_g", "bendoherb_g"
    ),
    labels = c(
      "b0_g" = "Grand Mean",
      "bendo_g" = "Endophyte",
      "bherb_g" = "Herbivory",
      "bsizesex_g" = "size:sex",
      "btemp_g" = "Temperature",
      "bendotemp_g" = "Endo:Temp",
      "bherbtemp_g" = "Herb:Temp",
      "bendoherb_g" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Growth)") +
  xlim(-2, 2) +
  ggtitle("B") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

posterior_flow <- as.array(fit_allsites_flow_temp_cv)
flow <- mcmc_intervals(posterior_flow, pars = quote_bare(
  b0_f, bendo_f, bherb_f, btemp_f,
  bendotemp_f, bherbtemp_f, bendoherb_f
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_f", "bendo_f", "bherb_f", "btemp_f",
      "bendotemp_f", "bherbtemp_f", "bendoherb_f"
    ),
    labels = c(
      "b0_f" = "Grand Mean",
      "bendo_f" = "Endophyte",
      "bherb_f" = "Herbivory",
      "bsizesex_f" = "size:sex",
      "btemp_f" = "Temperature",
      "bendotemp_f" = "Endo:Temp",
      "bherbtemp_f" = "Herb:Temp",
      "bendoherb_f" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Inflorescence)") +
  xlim(-6, 6) +
  ggtitle("C") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

posterior_flow <- as.array(fit_allsites_flow_temp_cv)

posterior_spi <- as.array(fit_allsites_spi_temp_cv)
spi <- mcmc_intervals(posterior_spi, pars = quote_bare(
  b0_spk, bendo_spk, bherb_spk, btemp_spk,
  bendotemp_spk, bherbtemp_spk, bendoherb_spk
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_spk", "bendo_spk", "bherb_spk", "btemp_spk",
      "bendotemp_spk", "bherbtemp_spk", "bendoherb_spk"
    ),
    labels = c(
      "b0_spk" = "Grand Mean",
      "bendo_spk" = "Endophyte",
      "bherb_spk" = "Herbivory",
      "bsizesex_spk" = "size:sex",
      "btemp_spk" = "Temperature",
      "bendotemp_spk" = "Endo:Temp",
      "bherbtemp_spk" = "Herb:Temp",
      "bendoherb_spk" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Spikelet)") +
  xlim(-5, 5) +
  ggtitle("D") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Posterior_mean.pdf", useDingbats = F, height = 7, width = 8)
ggarrange(surv, grow, flow, spi + rremove("ylab"), ncol = 2, nrow = 2)
dev.off()

## Temp mean
posterior_surv <- as.array(fit_allsites_surv_temp_mean)
# color_scheme_set("red")
surv <- mcmc_intervals(posterior_surv, pars = quote_bare(
  b0_s, bendo_s, bherb_s, btemp_s,
  bendotemp_s, bherbtemp_s, bendoherb_s
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_s", "bendo_s", "bherb_s", "btemp_s",
      "bendotemp_s", "bherbtemp_s", "bendoherb_s"
    ),
    labels = c(
      "b0_s" = "Grand Mean",
      "bendo_s" = "Endophyte",
      "bherb_s" = "Herbivory",
      "bsizesex_s" = "size:sex",
      "btemp_s" = "Temperature",
      "bendotemp_s" = "Endo:Temp",
      "bherbtemp_s" = "Herb:Temp",
      "bendoherb_s" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Survival)") +
  xlim(-8, 8) +
  ggtitle("A") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

posterior_grow <- as.array(fit_allsites_grow_temp_mean)
grow <- mcmc_intervals(posterior_grow, pars = quote_bare(
  b0_g, bendo_g, bherb_g, btemp_g,
  bendotemp_g, bherbtemp_g, bendoherb_g
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_g", "bendo_g", "bherb_g", "btemp_g",
      "bendotemp_g", "bherbtemp_g", "bendoherb_g"
    ),
    labels = c(
      "b0_g" = "Grand Mean",
      "bendo_g" = "Endophyte",
      "bherb_g" = "Herbivory",
      "bsizesex_g" = "size:sex",
      "btemp_g" = "Temperature",
      "bendotemp_g" = "Endo:Temp",
      "bherbtemp_g" = "Herb:Temp",
      "bendoherb_g" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Growth)") +
  xlim(-2, 2) +
  ggtitle("B") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

posterior_flow <- as.array(fit_allsites_flow_temp_mean)
flow <- mcmc_intervals(posterior_flow, pars = quote_bare(
  b0_f, bendo_f, bherb_f, btemp_f,
  bendotemp_f, bherbtemp_f, bendoherb_f
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_f", "bendo_f", "bherb_f", "btemp_f",
      "bendotemp_f", "bherbtemp_f", "bendoherb_f"
    ),
    labels = c(
      "b0_f" = "Grand Mean",
      "bendo_f" = "Endophyte",
      "bherb_f" = "Herbivory",
      "bsizesex_f" = "size:sex",
      "btemp_f" = "Temperature",
      "bendotemp_f" = "Endo:Temp",
      "bherbtemp_f" = "Herb:Temp",
      "bendoherb_f" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Inflorescence)") +
  xlim(-6, 6) +
  ggtitle("C") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

posterior_spi <- as.array(fit_allsites_spi_temp_mean)
spi <- mcmc_intervals(posterior_spi, pars = quote_bare(
  b0_spk, bendo_spk, bherb_spk, btemp_spk,
  bendotemp_spk, bherbtemp_spk, bendoherb_spk
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_spk", "bendo_spk", "bherb_spk", "btemp_spk",
      "bendotemp_spk", "bherbtemp_spk", "bendoherb_spk"
    ),
    labels = c(
      "b0_spk" = "Grand Mean",
      "bendo_spk" = "Endophyte",
      "bherb_spk" = "Herbivory",
      "bsizesex_spk" = "size:sex",
      "btemp_spk" = "Temperature",
      "bendotemp_spk" = "Endo:Temp",
      "bherbtemp_spk" = "Herb:Temp",
      "bendoherb_spk" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Spikelet)") +
  xlim(-5, 5) +
  ggtitle("D") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Posterior_mean_temp_mean.pdf", useDingbats = F, height = 7, width = 8)
ggarrange(surv, grow, flow, spi + rremove("ylab"), ncol = 2, nrow = 2)
dev.off()

## Water CV
posterior_surv <- as.array(fit_allsites_surv_water_cv)
# color_scheme_set("red")
surv <- mcmc_intervals(posterior_surv, pars = quote_bare(
  b0_s, bendo_s, bherb_s, btemp_s,
  bendotemp_s, bherbtemp_s, bendoherb_s
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_s", "bendo_s", "bherb_s", "btemp_s",
      "bendotemp_s", "bherbtemp_s", "bendoherb_s"
    ),
    labels = c(
      "b0_s" = "Grand Mean",
      "bendo_s" = "Endophyte",
      "bherb_s" = "Herbivory",
      "bsizesex_s" = "size:sex",
      "btemp_s" = "Temperature",
      "bendotemp_s" = "Endo:Temp",
      "bherbtemp_s" = "Herb:Temp",
      "bendoherb_s" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Survival)") +
  xlim(-8, 8) +
  ggtitle("A") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

posterior_grow <- as.array(fit_allsites_grow_water_cv)
grow <- mcmc_intervals(posterior_grow, pars = quote_bare(
  b0_g, bendo_g, bherb_g, btemp_g,
  bendotemp_g, bherbtemp_g, bendoherb_g
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_g", "bendo_g", "bherb_g", "btemp_g",
      "bendotemp_g", "bherbtemp_g", "bendoherb_g"
    ),
    labels = c(
      "b0_g" = "Grand Mean",
      "bendo_g" = "Endophyte",
      "bherb_g" = "Herbivory",
      "bsizesex_g" = "size:sex",
      "btemp_g" = "Temperature",
      "bendotemp_g" = "Endo:Temp",
      "bherbtemp_g" = "Herb:Temp",
      "bendoherb_g" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Growth)") +
  xlim(-2, 2) +
  ggtitle("B") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

posterior_flow <- as.array(fit_allsites_flow_water_cv)
flow <- mcmc_intervals(posterior_flow, pars = quote_bare(
  b0_f, bendo_f, bherb_f, btemp_f,
  bendotemp_f, bherbtemp_f, bendoherb_f
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_f", "bendo_f", "bherb_f", "btemp_f",
      "bendotemp_f", "bherbtemp_f", "bendoherb_f"
    ),
    labels = c(
      "b0_f" = "Grand Mean",
      "bendo_f" = "Endophyte",
      "bherb_f" = "Herbivory",
      "bsizesex_f" = "size:sex",
      "btemp_f" = "Temperature",
      "bendotemp_f" = "Endo:Temp",
      "bherbtemp_f" = "Herb:Temp",
      "bendoherb_f" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Inflorescence)") +
  xlim(-6, 6) +
  ggtitle("C") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

posterior_spi <- as.array(fit_allsites_spi_water_cv)
spi <- mcmc_intervals(posterior_spi, pars = quote_bare(
  b0_spk, bendo_spk, bherb_spk, btemp_spk,
  bendotemp_spk, bherbtemp_spk, bendoherb_spk
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_spk", "bendo_spk", "bherb_spk", "btemp_spk",
      "bendotemp_spk", "bherbtemp_spk", "bendoherb_spk"
    ),
    labels = c(
      "b0_spk" = "Grand Mean",
      "bendo_spk" = "Endophyte",
      "bherb_spk" = "Herbivory",
      "bsizesex_spk" = "size:sex",
      "btemp_spk" = "Temperature",
      "bendotemp_spk" = "Endo:Temp",
      "bherbtemp_spk" = "Herb:Temp",
      "bendoherb_spk" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Spikelet)") +
  xlim(-5, 5) +
  ggtitle("D") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Posterior_mean_water_cv.pdf", useDingbats = F, height = 7, width = 8)
ggarrange(surv, grow, flow, spi + rremove("ylab"), ncol = 2, nrow = 2)
dev.off()

## Water mean
posterior_surv <- as.array(fit_allsites_surv_water_mean)
# color_scheme_set("red")
surv <- mcmc_intervals(posterior_surv, pars = quote_bare(
  b0_s, bendo_s, bherb_s, btemp_s,
  bendotemp_s, bherbtemp_s, bendoherb_s
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_s", "bendo_s", "bherb_s", "btemp_s",
      "bendotemp_s", "bherbtemp_s", "bendoherb_s"
    ),
    labels = c(
      "b0_s" = "Grand Mean",
      "bendo_s" = "Endophyte",
      "bherb_s" = "Herbivory",
      "bsizesex_s" = "size:sex",
      "btemp_s" = "Temperature",
      "bendotemp_s" = "Endo:Temp",
      "bherbtemp_s" = "Herb:Temp",
      "bendoherb_s" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Survival)") +
  xlim(-8, 8) +
  ggtitle("A") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

posterior_grow <- as.array(fit_allsites_grow_water_mean)
grow <- mcmc_intervals(posterior_grow, pars = quote_bare(
  b0_g, bendo_g, bherb_g, btemp_g,
  bendotemp_g, bherbtemp_g, bendoherb_g
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_g", "bendo_g", "bherb_g", "btemp_g",
      "bendotemp_g", "bherbtemp_g", "bendoherb_g"
    ),
    labels = c(
      "b0_g" = "Grand Mean",
      "bendo_g" = "Endophyte",
      "bherb_g" = "Herbivory",
      "bsizesex_g" = "size:sex",
      "btemp_g" = "Temperature",
      "bendotemp_g" = "Endo:Temp",
      "bherbtemp_g" = "Herb:Temp",
      "bendoherb_g" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Growth)") +
  xlim(-2, 2) +
  ggtitle("B") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

posterior_flow <- as.array(fit_allsites_flow_water_mean)
flow <- mcmc_intervals(posterior_flow, pars = quote_bare(
  b0_f, bendo_f, bherb_f, btemp_f,
  bendotemp_f, bherbtemp_f, bendoherb_f
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_f", "bendo_f", "bherb_f", "btemp_f",
      "bendotemp_f", "bherbtemp_f", "bendoherb_f"
    ),
    labels = c(
      "b0_f" = "Grand Mean",
      "bendo_f" = "Endophyte",
      "bherb_f" = "Herbivory",
      "bsizesex_f" = "size:sex",
      "btemp_f" = "Temperature",
      "bendotemp_f" = "Endo:Temp",
      "bherbtemp_f" = "Herb:Temp",
      "bendoherb_f" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Inflorescence)") +
  xlim(-6, 6) +
  ggtitle("C") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

posterior_spi <- as.array(fit_allsites_spi_water_mean)
spi <- mcmc_intervals(posterior_spi, pars = quote_bare(
  b0_spk, bendo_spk, bherb_spk, btemp_spk,
  bendotemp_spk, bherbtemp_spk, bendoherb_spk
)) +
  ggplot2::scale_y_discrete(
    limits = c(
      "b0_spk", "bendo_spk", "bherb_spk", "btemp_spk",
      "bendotemp_spk", "bherbtemp_spk", "bendoherb_spk"
    ),
    labels = c(
      "b0_spk" = "Grand Mean",
      "bendo_spk" = "Endophyte",
      "bherb_spk" = "Herbivory",
      "bsizesex_spk" = "size:sex",
      "btemp_spk" = "Temperature",
      "bendotemp_spk" = "Endo:Temp",
      "bherbtemp_spk" = "Herb:Temp",
      "bendoherb_spk" = "Endo:Herb"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:") +
  xlab("Posterior estimates (Spikelet)") +
  xlim(-5, 5) +
  ggtitle("D") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr() +
  theme(
    axis.title.x = element_text(family = "Helvetica", colour = "black", size = 13),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.line.x = element_line(linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.5)
  )

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Posterior_mean_water_mean.pdf", useDingbats = F, height = 7, width = 8)
ggarrange(surv, grow, flow, spi + rremove("ylab"), ncol = 2, nrow = 2)
dev.off()

## re-format data for plotting
data_summary <- function(data, varname, groupnames) {
  require(plyr)
  summary_func <- function(x, col) {
    c(
      mean = mean(x[[col]], na.rm = TRUE),
      sd = sd(x[[col]], na.rm = TRUE)
    )
  }
  data_sum <- ddply(data, groupnames,
    .fun = summary_func,
    varname
  )
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
# Temperature
## Survival
demography_climate_elvi_surv %>%
  filter(Herbivory == "1") -> elvi_surv_no_herb
elvi_surv_no_herb_plot <- data_summary(elvi_surv_no_herb,
  varname = "surv_t1",
  groupnames = c("Endo", "temp_mean")
)
# Convert dose to a factor variable
elvi_surv_no_herb_plot$temp_mean <- round(elvi_surv_no_herb_plot$temp_mean, 2)
elvi_surv_no_herb_plot$temp_mean <- as.factor(elvi_surv_no_herb_plot$temp_mean)
elvi_surv_no_herb_plot$Endo <- as.factor(elvi_surv_no_herb_plot$Endo)

demography_climate_elvi_surv %>%
  filter(Herbivory == "2") -> elvi_surv_herb
elvi_surv_herb_plot <- data_summary(elvi_surv_herb,
  varname = "surv_t1",
  groupnames = c("Endo", "temp_mean")
)
# Convert dose to a factor variable
elvi_surv_herb_plot$temp_mean <- round(elvi_surv_herb_plot$temp_mean, 2)
elvi_surv_herb_plot$temp_mean <- as.factor(elvi_surv_herb_plot$temp_mean)
elvi_surv_herb_plot$Endo <- as.factor(elvi_surv_herb_plot$Endo)

sur_h <- ggplot(elvi_surv_no_herb_plot, aes(x = temp_mean, y = surv_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = surv_t1 - sd, ymax = surv_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Temperature (C)", y = "Survival") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = c(0.2, 0.2),
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )

surv <- ggplot(elvi_surv_herb_plot, aes(x = temp_mean, y = surv_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = surv_t1 - sd, ymax = surv_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Temperature (C)", y = "Growth") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )


## Growth
demography_climate_elvi_grow %>%
  filter(Herbivory == "1") -> elvi_grow_no_herb

demography_climate_elvi_grow %>%
  filter(Herbivory == "2") -> elvi_grow_herb

elvi_grow_no_herb_plot <- data_summary(elvi_grow_no_herb,
  varname = "grow",
  groupnames = c("Endo", "temp_mean")
)
# Convert dose to a factor variable
elvi_grow_no_herb_plot$temp_mean <- round(elvi_grow_no_herb_plot$temp_mean, 2)
elvi_grow_no_herb_plot$temp_mean <- as.factor(elvi_grow_no_herb_plot$temp_mean)
elvi_grow_no_herb_plot$Endo <- as.factor(elvi_grow_no_herb_plot$Endo)
elvi_grow_herb_plot <- data_summary(elvi_grow_herb,
  varname = "grow",
  groupnames = c("Endo", "temp_mean")
)
# Convert dose to a factor variable
elvi_grow_herb_plot$temp_mean <- round(elvi_grow_herb_plot$temp_mean, 2)
elvi_grow_herb_plot$temp_mean <- as.factor(elvi_grow_herb_plot$temp_mean)
elvi_grow_herb_plot$Endo <- as.factor(elvi_grow_herb_plot$Endo)
grow_h <- ggplot(elvi_grow_no_herb_plot, aes(x = temp_mean, y = grow, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = grow - sd, ymax = grow + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "Fenced", x = "Temperature (C)", y = "Growth") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = c(0.2, 0.23),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "grey90")
  )

grow <- ggplot(elvi_grow_herb_plot, aes(x = temp_mean, y = grow, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = grow - sd, ymax = grow + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "Unfenced", x = "Temperature (C)", y = "Growth") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )

## inflorescence
demography_climate_elvi_flowering %>%
  filter(Herbivory == "1") -> elvi_flow_no_herb
demography_climate_elvi_flowering %>%
  filter(Herbivory == "2") -> elvi_flow_herb

elvi_flow_no_herb_plot <- data_summary(elvi_flow_no_herb,
  varname = "flow_t1",
  groupnames = c("Endo", "temp_mean")
)
# Convert dose to a factor variable
elvi_flow_no_herb_plot$temp_mean <- round(elvi_flow_no_herb_plot$temp_mean, 2)
elvi_flow_no_herb_plot$temp_mean <- as.factor(elvi_flow_no_herb_plot$temp_mean)
elvi_flow_no_herb_plot$Endo <- as.factor(elvi_flow_no_herb_plot$Endo)
elvi_flow_herb_plot <- data_summary(elvi_flow_herb,
  varname = "flow_t1",
  groupnames = c("Endo", "temp_mean")
)
# Convert dose to a factor variable
elvi_flow_herb_plot$temp_mean <- round(elvi_flow_herb_plot$temp_mean, 2)
elvi_flow_herb_plot$temp_mean <- as.factor(elvi_flow_herb_plot$temp_mean)
elvi_flow_herb_plot$Endo <- as.factor(elvi_flow_herb_plot$Endo)

flow_h <- ggplot(elvi_flow_no_herb_plot, aes(x = temp_mean, y = flow_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = flow_t1 - sd, ymax = flow_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Temperature (C)", y = "# Inflorescences") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray95")
  )

flow <- ggplot(elvi_flow_herb_plot, aes(x = temp_mean, y = flow_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = flow_t1 - sd, ymax = flow_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Temperature (C)", y = " # Inflorescences") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )

## spikelet
demography_climate_elvi_spikelet %>%
  filter(Herbivory == "1") -> elvi_spi_no_herb
demography_climate_elvi_spikelet %>%
  filter(Herbivory == "2") -> elvi_spi_herb

elvi_spi_no_herb_plot <- data_summary(elvi_spi_no_herb,
  varname = "spikelet_t1",
  groupnames = c("Endo", "temp_mean")
)
# Convert dose to a factor variable
elvi_spi_no_herb_plot$temp_mean <- round(elvi_spi_no_herb_plot$temp_mean, 2)
elvi_spi_no_herb_plot$temp_mean <- as.factor(elvi_spi_no_herb_plot$temp_mean)
elvi_spi_no_herb_plot$Endo <- as.factor(elvi_spi_no_herb_plot$Endo)
elvi_spi_herb_plot <- data_summary(elvi_spi_herb,
  varname = "spikelet_t1",
  groupnames = c("Endo", "temp_mean")
)
# Convert dose to a factor variable
elvi_spi_herb_plot$temp_mean <- round(elvi_spi_herb_plot$temp_mean, 2)
elvi_spi_herb_plot$temp_mean <- as.factor(elvi_spi_herb_plot$temp_mean)
elvi_spi_herb_plot$Endo <- as.factor(elvi_spi_herb_plot$Endo)

spi_h <- ggplot(elvi_spi_no_herb_plot, aes(x = temp_mean, y = spikelet_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = spikelet_t1 - sd, ymax = spikelet_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Temperature (C)", y = "# Spikelets") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )

spi <- ggplot(elvi_spi_herb_plot, aes(x = temp_mean, y = spikelet_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = spikelet_t1 - sd, ymax = spikelet_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Temperature (C)", y = "# Spikelets") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Fig_temp_mean.pdf", useDingbats = F, height = 8.5, width = 6)
ggarrange(grow_h, grow + rremove("ylab"), flow_h, flow + rremove("ylab"), spi_h, spi + rremove("ylab"), ncol = 2, nrow = 3)
dev.off()


# Temperature CGV
## Survival
elvi_surv_no_herb_plot_t_cv <- data_summary(elvi_surv_no_herb,
  varname = "surv_t1",
  groupnames = c("Endo", "temp_cv")
)
# Convert dose to a factor variable
elvi_surv_no_herb_plot_t_cv$temp_cv <- round(elvi_surv_no_herb_plot_t_cv$temp_cv, 2)
elvi_surv_no_herb_plot_t_cv$temp_cv <- as.factor(elvi_surv_no_herb_plot_t_cv$temp_cv)
elvi_surv_no_herb_plot_t_cv$Endo <- as.factor(elvi_surv_no_herb_plot_t_cv$Endo)

elvi_surv_herb_plot_t_cv <- data_summary(elvi_surv_herb,
  varname = "surv_t1",
  groupnames = c("Endo", "temp_cv")
)
# Convert dose to a factor variable
elvi_surv_herb_plot_t_cv$temp_cv <- round(elvi_surv_herb_plot_t_cv$temp_cv, 2)
elvi_surv_herb_plot_t_cv$temp_cv <- as.factor(elvi_surv_herb_plot_t_cv$temp_cv)
elvi_surv_herb_plot_t_cv$Endo <- as.factor(elvi_surv_herb_plot_t_cv$Endo)

sur_h_t_cv <- ggplot(elvi_surv_no_herb_plot_t_cv, aes(x = temp_cv, y = surv_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = surv_t1 - sd, ymax = surv_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Temperature (CV)", y = "Survival") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = c(0.2, 0.2),
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )

surv_t_cv <- ggplot(elvi_surv_herb_plot_t_cv, aes(x = temp_cv, y = surv_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = surv_t1 - sd, ymax = surv_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Temperature (CV)", y = "Growth") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )


## Growth
elvi_grow_no_herb_plot_t_cv <- data_summary(elvi_grow_no_herb,
  varname = "grow",
  groupnames = c("Endo", "temp_cv")
)
# Convert dose to a factor variable
elvi_grow_no_herb_plot_t_cv$temp_cv <- round(elvi_grow_no_herb_plot_t_cv$temp_cv, 2)
elvi_grow_no_herb_plot_t_cv$temp_cv <- as.factor(elvi_grow_no_herb_plot_t_cv$temp_cv)
elvi_grow_no_herb_plot_t_cv$Endo <- as.factor(elvi_grow_no_herb_plot_t_cv$Endo)

elvi_grow_herb_plot_t_cv <- data_summary(elvi_grow_herb,
  varname = "grow",
  groupnames = c("Endo", "temp_cv")
)
# Convert dose to a factor variable
elvi_grow_herb_plot_t_cv$temp_cv <- round(elvi_grow_herb_plot_t_cv$temp_cv, 2)
elvi_grow_herb_plot_t_cv$temp_cv <- as.factor(elvi_grow_herb_plot_t_cv$temp_cv)
elvi_grow_herb_plot_t_cv$Endo <- as.factor(elvi_grow_herb_plot_t_cv$Endo)

grow_h_t_cv <- ggplot(elvi_grow_no_herb_plot_t_cv, aes(x = temp_cv, y = grow, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = grow - sd, ymax = grow + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "Fenced", x = "Temperature (CV)", y = "Growth") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = c(0.2, 0.23),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "grey90")
  )

grow_t_cv <- ggplot(elvi_grow_herb_plot_t_cv, aes(x = temp_cv, y = grow, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = grow - sd, ymax = grow + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "Unfenced", x = "Temperature (CV)", y = "Growth") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )

## inflorescence
elvi_flow_no_herb_plot_t_cv <- data_summary(elvi_flow_no_herb,
  varname = "flow_t1",
  groupnames = c("Endo", "temp_cv")
)
# Convert dose to a factor variable
elvi_flow_no_herb_plot_t_cv$temp_cv <- round(elvi_flow_no_herb_plot_t_cv$temp_cv, 2)
elvi_flow_no_herb_plot_t_cv$temp_cv <- as.factor(elvi_flow_no_herb_plot_t_cv$temp_cv)
elvi_flow_no_herb_plot_t_cv$Endo <- as.factor(elvi_flow_no_herb_plot_t_cv$Endo)

elvi_flow_herb_plot_t_cv <- data_summary(elvi_flow_herb,
  varname = "flow_t1",
  groupnames = c("Endo", "temp_cv")
)
# Convert dose to a factor variable
elvi_flow_herb_plot_t_cv$temp_cv <- round(elvi_flow_herb_plot_t_cv$temp_cv, 2)
elvi_flow_herb_plot_t_cv$temp_cv <- as.factor(elvi_flow_herb_plot_t_cv$temp_cv)
elvi_flow_herb_plot_t_cv$Endo <- as.factor(elvi_flow_herb_plot_t_cv$Endo)

flow_h_cv <- ggplot(elvi_flow_herb_plot_t_cv, aes(x = temp_cv, y = flow_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = flow_t1 - sd, ymax = flow_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Temperature (CV)", y = "# Inflorescences") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray95")
  )

flow_cv <- ggplot(elvi_flow_herb_plot_t_cv, aes(x = temp_cv, y = flow_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = flow_t1 - sd, ymax = flow_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Temperature (CV)", y = " # Inflorescences") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )

## spikelet

elvi_spi_no_herb_plot_t_cv <- data_summary(elvi_spi_no_herb,
  varname = "spikelet_t1",
  groupnames = c("Endo", "temp_cv")
)
# Convert dose to a factor variable
elvi_spi_no_herb_plot_t_cv$temp_cv <- round(elvi_spi_no_herb_plot_t_cv$temp_cv, 2)
elvi_spi_no_herb_plot_t_cv$temp_cv <- as.factor(elvi_spi_no_herb_plot_t_cv$temp_cv)
elvi_spi_no_herb_plot_t_cv$Endo <- as.factor(elvi_spi_no_herb_plot_t_cv$Endo)

elvi_spi_herb_plot_t_cv <- data_summary(elvi_spi_herb,
  varname = "spikelet_t1",
  groupnames = c("Endo", "temp_cv")
)
# Convert dose to a factor variable
elvi_spi_herb_plot_t_cv$temp_cv <- round(elvi_spi_herb_plot_t_cv$temp_cv, 2)
elvi_spi_herb_plot_t_cv$temp_cv <- as.factor(elvi_spi_herb_plot_t_cv$temp_cv)
elvi_spi_herb_plot_t_cv$Endo <- as.factor(elvi_spi_herb_plot_t_cv$Endo)

spi_h_t_cv <- ggplot(elvi_spi_no_herb_plot_t_cv, aes(x = temp_cv, y = spikelet_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = spikelet_t1 - sd, ymax = spikelet_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Temperature (CV)", y = "# Spikelets") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )

spi_t_cv <- ggplot(elvi_spi_herb_plot_t_cv, aes(x = temp_cv, y = spikelet_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = spikelet_t1 - sd, ymax = spikelet_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Temperature (CV)", y = "# Spikelets") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Fig_temp_cv.pdf", useDingbats = F, height = 8.5, width = 6)
ggarrange(grow_h_t_cv, grow_t_cv + rremove("ylab"), flow_h_cv, flow_cv + rremove("ylab"), spi_h_t_cv, spi_h_t_cv + rremove("ylab"), ncol = 2, nrow = 3)
dev.off()

# Soil moisture
## Survival
elvi_surv_no_herb_plot_w <- data_summary(elvi_surv_no_herb,
  varname = "surv_t1",
  groupnames = c("Endo", "water_mean")
)
# Convert dose to a factor variable
elvi_surv_no_herb_plot_w$water_mean <- round(elvi_surv_no_herb_plot_w$water_mean, 2)
elvi_surv_no_herb_plot_w$water_mean <- as.factor(elvi_surv_no_herb_plot_w$water_mean)
elvi_surv_no_herb_plot_w$Endo <- as.factor(elvi_surv_no_herb_plot_w$Endo)

elvi_surv_herb_plot_w <- data_summary(elvi_surv_herb,
  varname = "surv_t1",
  groupnames = c("Endo", "water_mean")
)
# Convert dose to a factor variable
elvi_surv_herb_plot_w$water_mean <- round(elvi_surv_herb_plot_w$water_mean, 2)
elvi_surv_herb_plot_w$water_mean <- as.factor(elvi_surv_herb_plot_w$water_mean)
elvi_surv_herb_plot_w$Endo <- as.factor(elvi_surv_herb_plot_w$Endo)

sur_h_w <- ggplot(elvi_surv_no_herb_plot_w, aes(x = water_mean, y = surv_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = surv_t1 - sd, ymax = surv_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Soil moisture (%) ", y = "Survival") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = c(0.2, 0.2),
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )

surv_w <- ggplot(elvi_surv_herb_plot_w, aes(x = water_mean, y = surv_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = surv_t1 - sd, ymax = surv_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Soil moisture (%)", y = "Growth") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )


## Growth
elvi_grow_no_herb_plot_w <- data_summary(elvi_grow_no_herb,
  varname = "grow",
  groupnames = c("Endo", "water_mean")
)
# Convert dose to a factor variable
elvi_grow_no_herb_plot_w$water_mean <- round(elvi_grow_no_herb_plot_w$water_mean, 2)
elvi_grow_no_herb_plot_w$water_mean <- as.factor(elvi_grow_no_herb_plot_w$water_mean)
elvi_grow_no_herb_plot_w$Endo <- as.factor(elvi_grow_no_herb_plot_w$Endo)

elvi_grow_herb_plot_w <- data_summary(elvi_grow_herb,
  varname = "grow",
  groupnames = c("Endo", "water_mean")
)
# Convert dose to a factor variable
elvi_grow_herb_plot_w$water_mean <- round(elvi_grow_herb_plot_w$water_mean, 2)
elvi_grow_herb_plot_w$water_mean <- as.factor(elvi_grow_herb_plot_w$water_mean)
elvi_grow_herb_plot_w$Endo <- as.factor(elvi_grow_herb_plot_w$Endo)
grow_h_w <- ggplot(elvi_grow_no_herb_plot_w, aes(x = water_mean, y = grow, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = grow - sd, ymax = grow + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "Fenced", x = "TSoil moisture (%)", y = "Growth") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = c(0.4, 0.23),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "grey90")
  )

grow_w <- ggplot(elvi_grow_herb_plot_w, aes(x = water_mean, y = grow, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = grow - sd, ymax = grow + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "Unfenced", x = "Soil moisture (%)", y = "Growth") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )

## inflorescence
elvi_flow_no_herb_plot_w <- data_summary(elvi_flow_no_herb,
  varname = "flow_t1",
  groupnames = c("Endo", "water_mean")
)
# Convert dose to a factor variable
elvi_flow_no_herb_plot_w$water_mean <- round(elvi_flow_no_herb_plot_w$water_mean, 2)
elvi_flow_no_herb_plot_w$water_mean <- as.factor(elvi_flow_no_herb_plot_w$water_mean)
elvi_flow_no_herb_plot_w$Endo <- as.factor(elvi_flow_no_herb_plot_w$Endo)

elvi_flow_herb_plot_w <- data_summary(elvi_flow_herb,
  varname = "flow_t1",
  groupnames = c("Endo", "water_mean")
)
# Convert dose to a factor variable
elvi_flow_herb_plot_w$water_mean <- round(elvi_flow_herb_plot_w$water_mean, 2)
elvi_flow_herb_plot_w$water_mean <- as.factor(elvi_flow_herb_plot_w$water_mean)
elvi_flow_herb_plot_w$Endo <- as.factor(elvi_flow_herb_plot_w$Endo)

flow_h_w <- ggplot(elvi_flow_no_herb_plot_w, aes(x = water_mean, y = flow_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = flow_t1 - sd, ymax = flow_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Soil moisture (%)", y = "# Inflorescences") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray95")
  )

flow_w <- ggplot(elvi_flow_herb_plot_w, aes(x = water_mean, y = flow_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = flow_t1 - sd, ymax = flow_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Soil moisture (%)", y = " # Inflorescences") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )

## spikelet
elvi_spi_no_herb_plot_w <- data_summary(elvi_spi_no_herb,
  varname = "spikelet_t1",
  groupnames = c("Endo", "water_mean")
)
# Convert dose to a factor variable
elvi_spi_no_herb_plot_w$water_mean <- round(elvi_spi_no_herb_plot_w$water_mean, 2)
elvi_spi_no_herb_plot_w$water_mean <- as.factor(elvi_spi_no_herb_plot_w$water_mean)
elvi_spi_no_herb_plot_w$Endo <- as.factor(elvi_spi_no_herb_plot_w$Endo)
elvi_spi_herb_plot_w <- data_summary(elvi_spi_herb,
  varname = "spikelet_t1",
  groupnames = c("Endo", "water_mean")
)
# Convert dose to a factor variable
elvi_spi_herb_plot_w$water_mean <- round(elvi_spi_herb_plot_w$water_mean, 2)
elvi_spi_herb_plot_w$water_mean <- as.factor(elvi_spi_herb_plot_w$water_mean)
elvi_spi_herb_plot_w$Endo <- as.factor(elvi_spi_herb_plot_w$Endo)

spi_h_w <- ggplot(elvi_spi_no_herb_plot_w, aes(x = water_mean, y = spikelet_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = spikelet_t1 - sd, ymax = spikelet_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Soil moisture (%)", y = "# Spikelets") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )

spi_w <- ggplot(elvi_spi_herb_plot_w, aes(x = water_mean, y = spikelet_t1, group = Endo, color = Endo)) +
  geom_line() +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin = spikelet_t1 - sd, ymax = spikelet_t1 + sd),
    width = .3,
    position = position_dodge(0.03)
  ) +
  labs(title = "", x = "Soil moisture (%)", y = "# Spikelets") +
  theme_bw() +
  scale_color_manual(name = "Endophyte", labels = c("1" = expression(E^"-"), "2" = expression(E^"+")), values = c("#999999", "#E69F00")) +
  theme(
    legend.position = "none",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", linetype = "blank", fill = "gray90")
  )

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Fig_water_mean.pdf", useDingbats = F, height = 8.5, width = 6)
ggarrange(grow_h_w, grow_w + rremove("ylab"), flow_h_w, flow_w + rremove("ylab"), spi_h_w, spi_w + rremove("ylab"), ncol = 2, nrow = 3)
dev.off()
