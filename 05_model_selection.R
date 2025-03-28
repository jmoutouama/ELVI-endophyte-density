# Project:
# Purpose: Model comparison based on epdl/looic 
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

# Survival----
fit_surv_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/hi11gxhpqlrdfg389ir0w/fit_surv_ppt.rds?rlkey=22ujyjnm74c6pw9biw50uasvu&dl=1"))
fit_surv_spei <- readRDS(url("https://www.dropbox.com/scl/fi/0js0md2myjvl2scu69bnm/fit_surv_spei.rds?rlkey=scn11z3a3epfgis8y8jrxke91&dl=1"))
fit_surv_distance <- readRDS(url("https://www.dropbox.com/scl/fi/gxc8edjzdjvsrtlb8zm5o/fit_surv_distance.rds?rlkey=bmtq9q0bxf6ooafq8pz3tyavp&dl=1"))

log_lik_surv_ppt <- loo::extract_log_lik(fit_surv_ppt, merge_chains = FALSE)
r_eff_surv_ppt <- loo::relative_eff(exp(log_lik_surv_ppt))
loo_surv_ppt <- loo(log_lik_surv_ppt, r_eff = r_eff_surv_ppt, cores = 4)
plot(loo_surv_ppt)

log_lik_surv_spei <- loo::extract_log_lik(fit_surv_spei, merge_chains = FALSE)
r_eff_surv_spei <- loo::relative_eff(exp(log_lik_surv_spei))
loo_surv_spei <- loo(log_lik_surv_spei, r_eff = r_eff_surv_spei, cores = 4)
plot(loo_surv_spei)

log_lik_surv_distance <- loo::extract_log_lik(fit_surv_distance, merge_chains = FALSE)
r_eff_surv_distance <- loo::relative_eff(exp(log_lik_surv_distance))
loo_surv_distance <- loo(log_lik_surv_distance, r_eff = r_eff_surv_distance, cores = 4)
plot(loo_surv_distance)

(comp_surv <- loo::loo_compare(loo_surv_ppt,loo_surv_spei, loo_surv_distance))

# Growth----
fit_grow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/85pzzrgvkxwogoybr004t/fit_grow_ppt.rds?rlkey=rz2xlg00u1aqhkxeaix7wiu9e&dl=1"))
fit_grow_spei <- readRDS(url("https://www.dropbox.com/scl/fi/4iuz5ay461qkjv732b5yb/fit_grow_spei.rds?rlkey=gq68ixyetb3v4o3ds7uo24cwk&dl=1"))
fit_grow_distance <- readRDS(url("https://www.dropbox.com/scl/fi/fhwjeizspvvd2dbz11195/fit_grow_distance.rds?rlkey=yt863sjawv1zliwem6a9ved1h&dl=1"))

log_lik_grow_ppt <- loo::extract_log_lik(fit_grow_ppt, merge_chains = FALSE)
r_eff_grow_ppt <- loo::relative_eff(exp(log_lik_grow_ppt))
loo_grow_ppt <- loo(log_lik_grow_ppt, r_eff = r_eff_grow_ppt, cores = 4)
plot(loo_grow_ppt)

log_lik_grow_spei <- loo::extract_log_lik(fit_grow_spei, merge_chains = FALSE)
r_eff_grow_spei <- loo::relative_eff(exp(log_lik_grow_spei))
loo_grow_spei <- loo(log_lik_grow_spei, r_eff = r_eff_grow_spei, cores = 4)
plot(loo_grow_spei)

log_lik_grow_distance <- loo::extract_log_lik(fit_grow_distance, merge_chains = FALSE)
r_eff_grow_distance <- loo::relative_eff(exp(log_lik_grow_distance))
loo_grow_distance <- loo(log_lik_grow_distance, r_eff = r_eff_grow_distance, cores = 4)
plot(loo_grow_distance)

(comp_grow <- loo::loo_compare(loo_grow_ppt,loo_grow_spei, loo_grow_distance))

# Flowering----
fit_flow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/ajo9euxdtxfmle2g006l3/fit_flow_ppt.rds?rlkey=78cknj1yaqs2m5cyicdzg8up4&dl=1"))
fit_flow_spei <- readRDS(url("https://www.dropbox.com/scl/fi/p4t60pr13u5qtt6hnz8kv/fit_flow_spei.rds?rlkey=5on6zasr6c1qt9wzaqux7lk0d&dl=1"))
fit_flow_distance <- readRDS(url("https://www.dropbox.com/scl/fi/y7qvz4pmy2t2j00gqvwqd/fit_flow_distance.rds?rlkey=qe47fg37dpi4z2lu1ur1lsbkh&dl=1"))

log_lik_flow_ppt <- loo::extract_log_lik(fit_flow_ppt, merge_chains = FALSE)
r_eff_flow_ppt <- loo::relative_eff(exp(log_lik_flow_ppt))
loo_flow_ppt <- loo(log_lik_flow_ppt, r_eff = r_eff_flow_ppt, cores = 4)
# plot(loo_flow_ppt)

log_lik_flow_spei <- loo::extract_log_lik(fit_flow_spei, merge_chains = FALSE)
r_eff_flow_spei <- loo::relative_eff(exp(log_lik_flow_spei))
loo_flow_spei <- loo(log_lik_flow_spei, r_eff = r_eff_flow_spei, cores = 4)
# plot(loo_flow_spei)

log_lik_flow_distance <- loo::extract_log_lik(fit_flow_distance, merge_chains = FALSE)
r_eff_flow_distance <- loo::relative_eff(exp(log_lik_flow_distance))
loo_flow_distance <- loo(log_lik_flow_distance, r_eff = r_eff_flow_distance, cores = 4)
# plot(loo_flow_distance)

(comp_flow <- loo::loo_compare(loo_flow_ppt,loo_flow_spei, loo_flow_distance))


# Spikelet----
fit_spik_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/9klrwl866stup2x3m2p9x/fit_spik_ppt.rds?rlkey=okno4gvn9dj5q8uhpnu36f59i&dl=1"))
fit_spik_spei <- readRDS(url("https://www.dropbox.com/scl/fi/uj8ok9u0hp2uw6gnpqhv3/fit_spik_spei.rds?rlkey=wfbjdxnjjx2qetc7cg8rinbg2&dl=1"))
fit_spik_distance <- readRDS(url("https://www.dropbox.com/scl/fi/1m01hqgokktxz77trnuu3/fit_spik_distance.rds?rlkey=jxx4nt8qgwgiospq870anu9mn&dl=1"))

log_lik_spi_ppt <- loo::extract_log_lik(fit_spik_ppt, merge_chains = FALSE)
r_eff_spi_ppt <- loo::relative_eff(exp(log_lik_spi_ppt))
loo_spi_ppt <- loo(log_lik_spi_ppt, r_eff = r_eff_spi_ppt, cores = 4)
#plot(loo_spi_ppt)

log_lik_spi_spei <- loo::extract_log_lik(fit_spik_spei, merge_chains = FALSE)
r_eff_spi_spei <- loo::relative_eff(exp(log_lik_spi_spei))
loo_spi_spei <- loo(log_lik_spi_spei, r_eff = r_eff_spi_spei, cores = 4)
plot(loo_spi_spei)

log_lik_spi_distance <- loo::extract_log_lik(fit_spik_distance, merge_chains = FALSE)
r_eff_spi_distance <- loo::relative_eff(exp(log_lik_spi_distance))
loo_spi_distance <- loo(log_lik_spi_distance, r_eff = r_eff_spi_distance, cores = 4)
plot(loo_spi_distance)


(comp_spi <- loo::loo_compare(loo_spi_temp_mean,loo_spi_temp_cv, loo_spi_water_mean,loo_spi_water_cv))
