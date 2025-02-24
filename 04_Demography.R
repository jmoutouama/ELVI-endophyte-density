# Project: 
# Purpose: Fit vital rate models to test the effect of grass-endophyte symbiosis and endophyte hyphal density on  vital rate models (survival, growth, flowering,fertility).
# Note: Raster files are too large to provide in public repository. They are stored on a local machine
# Authors: Jacob Moutouama
# Date last modified (Y-M-D): 2024-08-03
rm(list = ls())
# load packages
#remove.packages(c("StanHeaders", "rstan"))
#install.packages("rstan", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
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
#library(countreg)
library(rmutil)
library(actuar)
#library(SPEI)
library(LaplacesDemon)
library(ggpubr)
library(raster)
#library(rgdal)
library(readxl)
library(ggsci)
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("scater")
# library(scater)
library(BiocManager)
library(swfscMisc)

# Define some basic functions that we'll use later 
quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply(deparse )
}

invlogit<-function(x){exp(x)/(1+exp(x))}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# Demographic data -----
# Merge the demographic census
datini<-read.csv("https://www.dropbox.com/scl/fi/exwmw8z8vp1qkf8inyeoq/Initialdata.csv?rlkey=kez08s92dgh9v0i08269kx1iq&dl=1", stringsAsFactors = F)
dat23<-read.csv("https://www.dropbox.com/scl/fi/9ob0vpu2xdq8x7u48866s/census2023.csv?rlkey=i2loj3fezymq1p41bo5lsj3uj&dl=1", stringsAsFactors = F)
dat24<-read.csv("https://www.dropbox.com/scl/fi/s8pnf1j7c85g6jwc944vw/census2024.csv?rlkey=kwt2x8k16q4w7gndm42komj6o&dl=1", stringsAsFactors = F)
datherbivory<-read.csv("https://www.dropbox.com/scl/fi/suy4twdhy36el0k7ytqsi/herbivory.csv?rlkey=hs4xbjn1zrpnpitry30ng538d&dl=1", stringsAsFactors = F)
# unique(datini$Site)
# unique(datini$dat23)
# unique(datini$dat24)
# names(dat23)
# calculate the total spikelet for each census
dat23 %>% 
  mutate(spikelet_23=round(rowMeans(across(Spikelet_A:Spikelet_C),na.rm=T)),digit=0)->dat23_spike
dat24 %>% 
  mutate(spikelet_24=round(rowMeans(across(Spikelet_A:Spikelet_C),na.rm=T),digit=0),Inf_24=round(rowMeans(across(attachedInf_24:brokenInf_24),na.rm=T),digit=0))->dat24_spike

## Merge the initial data with the 23 data and the 23 data with the 24 -----
datini23 <- left_join(x = datini,y =dat23_spike,by=c("Tag_ID"))
names(datini23)
dat2324 <- left_join(x = datini23 ,y =dat24_spike,by=c("Tag_ID")) 
names(dat2324)
dat2324 %>% 
  mutate(tiller_t=Tiller_23,
         tiller_t1=Tiller_24,
         inf_t=Inf_23,
         inf_t1=Inf_24,
         spikelet_t=spikelet_23,
         spikelet_t1=spikelet_24,
         tiller_Herb_t=tiller_Herb,
         tiller_Herb_t1=tiller_herb_24) %>% 
  dplyr::select(Site,
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
                date_24)->dat2324_t_t1
#names(dat2324_t_t1)
## Merge the demographic data with the herbivory data -----
dat2324_t_t1_herb<-left_join(x=dat2324_t_t1,y=datherbivory,by=c("Site","Plot","Species"))# Merge the demographic data with the herbivory data
head(dat2324_t_t1_herb)
unique(dat2324_t_t1_herb$Species)
#view(dat2324_t_t1_herb)

## Merge the demographic data with the climatic data -----
demography_climate<-left_join(x=dat2324_t_t1_herb_ELVI_clean,y=HOBO_summary_clean,by=c("Site"))# Merge the demographic data with the temperature data
demography_climate$surv1<-1*(!is.na(demography_climate$tiller_t) & !is.na(demography_climate$tiller_t1))
demography_climate$site_plot<-interaction(demography_climate$Site,demography_climate$Plot)
demography_climate$grow<-(log(demography_climate$tiller_t1+1) - log(demography_climate$tiller_t+1))# Relative growth rate
# names(demography_climate)
# view(demography_climate)
#summary(demography_climate)
hist(demography_climate$grow,main="")

# H1: We hypothesized that stress associated with aridity and low precipitation would strengthen the plant-fungal mutualism, such that the fitness benefits of endophyte symbiosis are maximized at the range edge. ----

## Survival----
## Read and format survival data to build the model
demography_climate %>% 
  subset( tiller_t > 0 )%>%
  dplyr::select(Population, Site, Plot,site_plot, Endo, Herbivory,
                tiller_t, surv1,temp_mean,temp_cv,water_mean,water_cv)%>% 
  na.omit %>% 
  mutate( Site= Site %>% as.factor %>% as.numeric,
          Plot = Plot %>% as.factor %>% as.numeric,
          site_plot=site_plot %>% as.factor %>% as.numeric,
          Endo = Endo %>% as.factor %>% as.numeric,
          Herbivory=Herbivory %>% as.factor %>% as.numeric,
          Population = Population %>% as.factor %>% as.numeric ) %>%
  mutate( log_size_t0 = log(tiller_t),
          surv_t1=surv1,
          log_temp_mean = log(temp_mean),
          log_temp_cv = log(temp_cv),
          log_water_mean = log(water_mean),
          log_water_cv = log(water_cv))->demography_climate_elvi_surv

## Separate each variable to use the same model stan
data_sites_surv_temp_mean <- list( n_sites    = demography_climate_elvi_surv$Site %>% n_distinct,
                           n_pops  = demography_climate_elvi_surv$Population %>% n_distinct(),
                           # survival data
                           n_plot_s = demography_climate_elvi_surv$Plot %>% n_distinct,
                           site_s   = demography_climate_elvi_surv$Site,
                           pop_s =  demography_climate_elvi_surv$Population,
                           plot_s  = demography_climate_elvi_surv$Plot,
                           temp_s=as.vector(demography_climate_elvi_surv$log_temp_mean),
                          endo_s  = demography_climate_elvi_surv$Endo-1,
                          herb_s  = demography_climate_elvi_surv$Herbivory-1,
                           size_s   = demography_climate_elvi_surv$log_size_t0,
                           y_s      = demography_climate_elvi_surv$surv_t1,
                           n_s      = nrow(demography_climate_elvi_surv))
data_sites_surv_temp_cv <- list( n_sites    = demography_climate_elvi_surv$Site %>% n_distinct,
                                   n_pops  = demography_climate_elvi_surv$Population %>% n_distinct(),
                                   # survival data
                                   n_plot_s = demography_climate_elvi_surv$Plot %>% n_distinct,
                                   site_s   = demography_climate_elvi_surv$Site,
                                   pop_s =  demography_climate_elvi_surv$Population,
                                   plot_s  = demography_climate_elvi_surv$Plot,
                                   temp_s=as.vector(demography_climate_elvi_surv$log_temp_cv),
                                   endo_s  = demography_climate_elvi_surv$Endo-1,
                                   herb_s  = demography_climate_elvi_surv$Herbivory-1,
                                   size_s   = demography_climate_elvi_surv$log_size_t0,
                                   y_s      = demography_climate_elvi_surv$surv_t1,
                                   n_s      = nrow(demography_climate_elvi_surv))
data_sites_surv_water_mean <- list( n_sites    = demography_climate_elvi_surv$Site %>% n_distinct,
                                 n_pops  = demography_climate_elvi_surv$Population %>% n_distinct(),
                                 # survival data
                                 n_plot_s = demography_climate_elvi_surv$Plot %>% n_distinct,
                                 site_s   = demography_climate_elvi_surv$Site,
                                 pop_s =  demography_climate_elvi_surv$Population,
                                 plot_s  = demography_climate_elvi_surv$Plot,
                                 temp_s=as.vector(demography_climate_elvi_surv$water_mean),
                                 endo_s  = demography_climate_elvi_surv$Endo-1,
                                 herb_s  = demography_climate_elvi_surv$Herbivory-1,
                                 size_s   = demography_climate_elvi_surv$log_size_t0,
                                 y_s      = demography_climate_elvi_surv$surv_t1,
                                 n_s      = nrow(demography_climate_elvi_surv))
data_sites_surv_water_cv <- list( n_sites    = demography_climate_elvi_surv$Site %>% n_distinct,
                                    n_pops  = demography_climate_elvi_surv$Population %>% n_distinct(),
                                    # survival data
                                    n_plot_s = demography_climate_elvi_surv$Plot %>% n_distinct,
                                    site_s   = demography_climate_elvi_surv$Site,
                                    pop_s =  demography_climate_elvi_surv$Population,
                                    plot_s  = demography_climate_elvi_surv$Plot,
                                    temp_s=as.vector(demography_climate_elvi_surv$water_cv),
                                    endo_s  = demography_climate_elvi_surv$Endo-1,
                                    herb_s  = demography_climate_elvi_surv$Herbivory-1,
                                    size_s   = demography_climate_elvi_surv$log_size_t0,
                                    y_s      = demography_climate_elvi_surv$surv_t1,
                                    n_s      = nrow(demography_climate_elvi_surv))
## Running the stan model

# sim_pars <- list(
#   warmup = 2000, 
#   iter = 8000, 
#   thin = 2, 
#   chains = 4
# )

# fit_allsites_surv_temp_mean <- stan(
#  file = "D:/stan/elvi_survival.stan",
#  data = data_sites_surv_temp_mean,
#  warmup = sim_pars$warmup,
#  iter = sim_pars$iter,
#  thin = sim_pars$thin,
#  chains = sim_pars$chains)
# 
# fit_allsites_surv_temp_cv <- stan(
#   file = "D:/stan/elvi_survival.stan",
#   data = data_sites_surv_temp_cv,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains)
# 
# fit_allsites_surv_water_mean <- stan(
#   file = "D:/stan/elvi_survival.stan",
#   data = data_sites_surv_water_mean,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains)
# 
# fit_allsites_surv_water_cv <- stan(
#   file = "D:/stan/elvi_survival.stan",
#   data = data_sites_surv_water_cv,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains)

## Save RDS file for further use
# saveRDS(fit_allsites_surv_temp_mean, 'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_surv_temp_mean.rds')
# saveRDS(fit_allsites_surv_temp_cv, 'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_surv_temp_cv.rds')
# saveRDS(fit_allsites_surv_water_mean, 'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_surv_water_mean.rds')
# saveRDS(fit_allsites_surv_water_cv, 'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_surv_water_cv.rds')

## Flowering----
demography_climate %>% 
  subset( tiller_t1 > 0 )%>%
  dplyr::select( Population, Site, Plot,site_plot, Endo, Herbivory,
                tiller_t1, inf_t1,temp_mean,temp_cv,water_mean,water_cv)%>% 
  na.omit %>% 
  mutate( Site= Site %>% as.factor %>% as.numeric,
          Plot = Plot %>% as.factor %>% as.numeric,
          site_plot=site_plot %>% as.factor %>% as.numeric,
          Endo = Endo %>% as.factor %>% as.numeric,
          Herbivory=Herbivory %>% as.factor %>% as.numeric,
          Population = Population %>% as.factor %>% as.numeric ) %>%
  mutate( flow_t1=inf_t1,
          log_temp_mean = log(temp_mean),
          log_temp_cv = log(temp_cv),
          log_water_mean = log(water_mean),
          log_water_cv = log(water_cv))->demography_climate_elvi_flowering
summary(demography_climate_elvi_flowering)
# hist(demography_climate_elvi_flowering$flow_t1)
## Separate each variable to use the same model stan
data_sites_flow_temp_mean <- list( n_sites    = demography_climate_elvi_flowering$Site %>% n_distinct,
                                   n_pops  = demography_climate_elvi_flowering$Population %>% n_distinct(),
                                   # survival data
                                   n_plot_f = demography_climate_elvi_flowering$Plot %>% n_distinct,
                                   site_f   = demography_climate_elvi_flowering$Site,
                                   pop_f =  demography_climate_elvi_flowering$Population,
                                   plot_f  = demography_climate_elvi_flowering$Plot,
                                   temp_f=as.vector(demography_climate_elvi_flowering$log_temp_mean),
                                   endo_f  = demography_climate_elvi_flowering$Endo-1,
                                   herb_f  = demography_climate_elvi_flowering$Herbivory-1,
                                   y_f      = demography_climate_elvi_flowering$flow_t1,
                                   n_f      = nrow(demography_climate_elvi_flowering))

data_sites_flow_temp_cv <- list( n_sites    = demography_climate_elvi_flowering$Site %>% n_distinct,
                                   n_pops  = demography_climate_elvi_flowering$Population %>% n_distinct(),
                                   # survival data
                                   n_plot_f = demography_climate_elvi_flowering$Plot %>% n_distinct,
                                   site_f   = demography_climate_elvi_flowering$Site,
                                   pop_f =  demography_climate_elvi_flowering$Population,
                                   plot_f  = demography_climate_elvi_flowering$Plot,
                                   temp_f=as.vector(demography_climate_elvi_flowering$log_temp_cv),
                                   endo_f  = demography_climate_elvi_flowering$Endo-1,
                                   herb_f  = demography_climate_elvi_flowering$Herbivory-1,
                                   y_f      = demography_climate_elvi_flowering$flow_t1,
                                   n_f      = nrow(demography_climate_elvi_flowering))

data_sites_flow_water_mean <- list( n_sites    = demography_climate_elvi_flowering$Site %>% n_distinct,
                                   n_pops  = demography_climate_elvi_flowering$Population %>% n_distinct(),
                                   # survival data
                                   n_plot_f = demography_climate_elvi_flowering$Plot %>% n_distinct,
                                   site_f   = demography_climate_elvi_flowering$Site,
                                   pop_f =  demography_climate_elvi_flowering$Population,
                                   plot_f  = demography_climate_elvi_flowering$Plot,
                                   temp_f=as.vector(demography_climate_elvi_flowering$log_water_mean),
                                   endo_f  = demography_climate_elvi_flowering$Endo-1,
                                   herb_f  = demography_climate_elvi_flowering$Herbivory-1,
                                   y_f      = demography_climate_elvi_flowering$flow_t1,
                                   n_f      = nrow(demography_climate_elvi_flowering))

data_sites_flow_water_cv <- list( n_sites    = demography_climate_elvi_flowering$Site %>% n_distinct,
                                 n_pops  = demography_climate_elvi_flowering$Population %>% n_distinct(),
                                 # survival data
                                 n_plot_f = demography_climate_elvi_flowering$Plot %>% n_distinct,
                                 site_f   = demography_climate_elvi_flowering$Site,
                                 pop_f =  demography_climate_elvi_flowering$Population,
                                 plot_f  = demography_climate_elvi_flowering$Plot,
                                 temp_f=as.vector(demography_climate_elvi_flowering$log_water_cv),
                                 endo_f  = demography_climate_elvi_flowering$Endo-1,
                                 herb_f  = demography_climate_elvi_flowering$Herbivory-1,
                                 y_f      = demography_climate_elvi_flowering$flow_t1,
                                 n_f      = nrow(demography_climate_elvi_flowering))

# fit_allsites_flow_temp_mean <- stan(
#   file = "D:/stan/elvi_flowering.stan",
#   data = data_sites_flow_temp_mean,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains)
# 
# fit_allsites_flow_temp_cv <- stan(
#   file = "D:/stan/elvi_flowering.stan",
#   data = data_sites_flow_temp_cv,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains)
# 
# fit_allsites_flow_water_mean <- stan(
#   file = "D:/stan/elvi_flowering.stan",
#   data = data_sites_flow_water_mean,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains)
# 
# fit_allsites_flow_water_cv <- stan(
#   file = "D:/stan/elvi_flowering.stan",
#   data = data_sites_flow_water_cv,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains)

# saveRDS(fit_allsites_flow_temp_mean, 'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_flow_temp_mean.rds')
# saveRDS(fit_allsites_flow_temp_cv, 'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_flow_temp_cv.rds')
# saveRDS(fit_allsites_flow_water_mean, 'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_flow_water_mean.rds')
# saveRDS(fit_allsites_flow_water_cv, 'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_flow_water_cv.rds')

## Spikelet----
demography_climate %>% 
  subset( spikelet_t1 > 0 & tiller_t1 > 0 )%>%
  dplyr::select(Population, Site, Plot,site_plot, Endo, Herbivory,
                tiller_t1, spikelet_t1,temp_mean,temp_cv,water_mean,water_cv)%>% 
  na.omit %>% 
  mutate( Site= Site %>% as.factor %>% as.numeric,
          Plot = Plot %>% as.factor %>% as.numeric,
          site_plot=site_plot %>% as.factor %>% as.numeric,
          Endo = Endo %>% as.factor %>% as.numeric,
          Herbivory=Herbivory %>% as.factor %>% as.numeric,
          Population = Population %>% as.factor %>% as.numeric ) %>%
  mutate( spi_t1=spikelet_t1,
          log_temp_mean = log(temp_mean),
          log_temp_cv = log(temp_cv),
          log_water_mean = log(water_mean),
          log_water_cv = log(water_cv))->demography_climate_elvi_spikelet


data_sites_spi_temp_mean <- list( n_sites    = demography_climate_elvi_spikelet$Site %>% n_distinct,
                                   n_pops  = demography_climate_elvi_spikelet$Population %>% n_distinct(),
                                   # survival data
                                   n_plot_spk = demography_climate_elvi_spikelet$Plot %>% n_distinct,
                                   site_spk   = demography_climate_elvi_spikelet$Site,
                                   pop_spk =  demography_climate_elvi_spikelet$Population,
                                   plot_spk  = demography_climate_elvi_spikelet$Plot,
                                   temp_spk=as.vector(demography_climate_elvi_spikelet$log_temp_mean),
                                   endo_spk  = demography_climate_elvi_spikelet$Endo-1,
                                   herb_spk  = demography_climate_elvi_spikelet$Herbivory-1,
                                   y_spk      = demography_climate_elvi_spikelet$spi_t1,
                                   n_spk      = nrow(demography_climate_elvi_spikelet))

data_sites_spi_temp_cv <- list( n_sites    = demography_climate_elvi_spikelet$Site %>% n_distinct,
                                 n_pops  = demography_climate_elvi_spikelet$Population %>% n_distinct(),
                                 # survival data
                                 n_plot_spk = demography_climate_elvi_spikelet$Plot %>% n_distinct,
                                 site_spk   = demography_climate_elvi_spikelet$Site,
                                 pop_spk =  demography_climate_elvi_spikelet$Population,
                                 plot_spk  = demography_climate_elvi_spikelet$Plot,
                                 temp_spk=as.vector(demography_climate_elvi_spikelet$log_temp_cv),
                                 endo_spk  = demography_climate_elvi_spikelet$Endo-1,
                                 herb_spk  = demography_climate_elvi_spikelet$Herbivory-1,
                                 y_spk      = demography_climate_elvi_spikelet$spi_t1,
                                 n_spk      = nrow(demography_climate_elvi_spikelet))

data_sites_spi_water_mean <- list( n_sites    = demography_climate_elvi_spikelet$Site %>% n_distinct,
                                    n_pops  = demography_climate_elvi_spikelet$Population %>% n_distinct(),
                                    # survival data
                                    n_plot_spk = demography_climate_elvi_spikelet$Plot %>% n_distinct,
                                    site_spk   = demography_climate_elvi_spikelet$Site,
                                    pop_spk =  demography_climate_elvi_spikelet$Population,
                                    plot_spk  = demography_climate_elvi_spikelet$Plot,
                                    temp_spk=as.vector(demography_climate_elvi_spikelet$log_water_mean),
                                    endo_spk  = demography_climate_elvi_spikelet$Endo-1,
                                    herb_spk  = demography_climate_elvi_spikelet$Herbivory-1,
                                    y_spk      = demography_climate_elvi_spikelet$spi_t1,
                                    n_spk      = nrow(demography_climate_elvi_spikelet))

data_sites_spi_water_cv <- list( n_sites    = demography_climate_elvi_spikelet$Site %>% n_distinct,
                                  n_pops  = demography_climate_elvi_spikelet$Population %>% n_distinct(),
                                  # survival data
                                  n_plot_spk = demography_climate_elvi_spikelet$Plot %>% n_distinct,
                                  site_spk   = demography_climate_elvi_spikelet$Site,
                                  pop_spk =  demography_climate_elvi_spikelet$Population,
                                  plot_spk  = demography_climate_elvi_spikelet$Plot,
                                  temp_spk=as.vector(demography_climate_elvi_spikelet$log_water_cv),
                                  endo_spk  = demography_climate_elvi_spikelet$Endo-1,
                                  herb_spk  = demography_climate_elvi_spikelet$Herbivory-1,
                                  y_spk      = demography_climate_elvi_spikelet$spi_t1,
                                  n_spk      = nrow(demography_climate_elvi_spikelet))

# fit_allsites_spi_temp_mean <- stan(
#   file = "D:/stan/elvi_spilkelet.stan",
#   data = data_sites_spi_temp_mean,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains)
# 
# fit_allsites_spi_temp_cv <- stan(
#   file = "D:/stan/elvi_spilkelet.stan",
#   data = data_sites_spi_temp_cv,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains)
# 
# fit_allsites_spi_water_mean <- stan(
#   file = "D:/stan/elvi_spilkelet.stan",
#   data = data_sites_spi_water_mean,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains)
# 
# fit_allsites_spi_water_cv <- stan(
#   file = "D:/stan/elvi_spilkelet.stan",
#   data = data_sites_spi_water_cv,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains)

#Save RDS file for further use
# saveRDS(fit_allsites_spi_temp_mean,'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_spi_temp_mean.rds')
# saveRDS(fit_allsites_spi_temp_cv, 'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_spi_temp_cv.rds')
# saveRDS(fit_allsites_spi_water_mean, 'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_spi_water_mean.rds')
# saveRDS(fit_allsites_spi_water_cv, 'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_spi_water_cv.rds')

## Growth----
## Read and format the growth data to build the model
#hist(demography_climate_elvi$grow)
demography_climate %>% 
  subset( tiller_t > 0 & tiller_t1 > 0)%>%
  dplyr::select(Population, Site, Plot,site_plot, Endo, Herbivory,
                tiller_t, grow,temp_mean,temp_cv,water_mean,water_cv)%>% 
  na.omit %>% 
  mutate( Site= Site %>% as.factor %>% as.numeric,
          Plot = Plot %>% as.factor %>% as.numeric,
          site_plot=site_plot %>% as.factor %>% as.numeric,
          Endo = Endo %>% as.factor %>% as.numeric,
          Herbivory=Herbivory %>% as.factor %>% as.numeric,
          Population = Population %>% as.factor %>% as.numeric ) %>%
  mutate( log_size_t0 = log(tiller_t),
          grow_t1=grow,
          log_temp_mean = log(temp_mean),
          log_temp_cv = log(temp_cv),
          log_water_mean = log(water_mean),
          log_water_cv = log(water_cv))->demography_climate_elvi_grow

## Separate each variable to use the same model stan
data_sites_grow_temp_mean <- list( n_sites    = demography_climate_elvi_grow$Site %>% n_distinct,
                                   n_pops  = demography_climate_elvi_grow$Population %>% n_distinct(),
                                   # survival data
                                   n_plot_g = demography_climate_elvi_grow$Plot %>% n_distinct,
                                   site_g   = demography_climate_elvi_grow$Site,
                                   pop_g =  demography_climate_elvi_grow$Population,
                                   plot_g  = demography_climate_elvi_grow$Plot,
                                   temp_g=as.vector(demography_climate_elvi_grow$log_temp_mean),
                                   endo_g  = demography_climate_elvi_grow$Endo-1,
                                   herb_g  = demography_climate_elvi_grow$Herbivory-1,
                                   size_g   = demography_climate_elvi_grow$log_size_t0,
                                   y_g      = demography_climate_elvi_grow$grow_t1,
                                   n_g      = nrow(demography_climate_elvi_grow))
# summary(data_sites_grow_temp_mean$y_g)
data_sites_grow_temp_cv <- list( n_sites    = demography_climate_elvi_grow$Site %>% n_distinct,
                                 n_pops  = demography_climate_elvi_grow$Population %>% n_distinct(),
                                 # survival data
                                 n_plot_g = demography_climate_elvi_grow$Plot %>% n_distinct,
                                 site_g   = demography_climate_elvi_grow$Site,
                                 pop_g =  demography_climate_elvi_grow$Population,
                                 plot_g  = demography_climate_elvi_grow$Plot,
                                 temp_g=as.vector(demography_climate_elvi_grow$log_temp_cv),
                                 endo_g  = demography_climate_elvi_grow$Endo-1,
                                 herb_g  = demography_climate_elvi_grow$Herbivory-1,
                                 size_g   = demography_climate_elvi_grow$log_size_t0,
                                 y_g      = demography_climate_elvi_grow$grow_t1,
                                 n_g      = nrow(demography_climate_elvi_grow))
data_sites_grow_water_mean <- list( n_sites    = demography_climate_elvi_grow$Site %>% n_distinct,
                                    n_pops  = demography_climate_elvi_grow$Population %>% n_distinct(),
                                    # survival data
                                    n_plot_g = demography_climate_elvi_grow$Plot %>% n_distinct,
                                    site_g   = demography_climate_elvi_grow$Site,
                                    pop_g =  demography_climate_elvi_grow$Population,
                                    plot_g  = demography_climate_elvi_grow$Plot,
                                    temp_g=as.vector(demography_climate_elvi_grow$water_mean),
                                    endo_g  = demography_climate_elvi_grow$Endo-1,
                                    herb_g  = demography_climate_elvi_grow$Herbivory-1,
                                    size_g   = demography_climate_elvi_grow$log_size_t0,
                                    y_g      = demography_climate_elvi_grow$grow_t1,
                                    n_g      = nrow(demography_climate_elvi_grow))
data_sites_grow_water_cv <- list( n_sites    = demography_climate_elvi_grow$Site %>% n_distinct,
                                  n_pops  = demography_climate_elvi_grow$Population %>% n_distinct(),
                                  # survival data
                                  n_plot_g = demography_climate_elvi_grow$Plot %>% n_distinct,
                                  site_g   = demography_climate_elvi_grow$Site,
                                  pop_g =  demography_climate_elvi_grow$Population,
                                  plot_g  = demography_climate_elvi_grow$Plot,
                                  temp_g=as.vector(demography_climate_elvi_grow$water_cv),
                                  endo_g  = demography_climate_elvi_grow$Endo-1,
                                  herb_g  = demography_climate_elvi_grow$Herbivory-1,
                                  size_g   = demography_climate_elvi_grow$log_size_t0,
                                  y_g      = demography_climate_elvi_grow$grow_t1,
                                  n_g      = nrow(demography_climate_elvi_grow))


# fit_allsites_grow_temp_mean <- stan(
#   file = "D:/stan/elvi_growth.stan",
#   data = data_sites_grow_temp_mean,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains)
# 
# fit_allsites_grow_temp_cv <- stan(
#   file = "D:/stan/elvi_growth.stan",
#   data = data_sites_grow_temp_cv,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains)
# 
# fit_allsites_grow_water_mean <- stan(
#   file = "D:/stan/elvi_growth.stan",
#   data = data_sites_grow_water_mean,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains)
# 
# fit_allsites_grow_water_cv <- stan(
#   file = "D:/stan/elvi_growth.stan",
#   data = data_sites_grow_water_cv,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains)

# Save RDS file for further use
# saveRDS(fit_allsites_grow_temp_mean, 'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_grow_temp_mean.rds')
# saveRDS(fit_allsites_grow_temp_cv, 'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_grow_temp_cv.rds')
# saveRDS(fit_allsites_grow_water_mean, 'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_grow_water_mean.rds')
# saveRDS(fit_allsites_grow_water_cv, 'C:/Users/jm200/Documents/ELVI Stan Output/fit_allsites_grow_water_cv.rds')

## load stan output 
fit_allsites_surv_temp_mean <- readRDS(url("https://www.dropbox.com/scl/fi/lmervp4u3uil63izyf4d7/fit_allsites_surv_temp_mean.rds?rlkey=8xljiw4tbefq902jtsnvhrbdp&dl=1"))
fit_allsites_surv_temp_cv <- readRDS(url("https://www.dropbox.com/scl/fi/f5s80kit3t3zxcj40mzk0/fit_allsites_surv_temp_cv.rds?rlkey=vvtgc5m54qdj59bc56m0cm2gu&dl=1"))
fit_allsites_surv_water_mean <- readRDS(url("https://www.dropbox.com/scl/fi/ag4mbqpv5dquy6b55brf2/fit_allsites_surv_water_mean.rds?rlkey=7kubpn2n1rc3zwiq40olrfrsf&dl=1"))
fit_allsites_surv_water_cv <- readRDS(url("https://www.dropbox.com/scl/fi/36xe54ctx1oiup8u1rqr6/fit_allsites_surv_water_cv.rds?rlkey=u0p31mhd7s959rxd9vwkn0zhw&dl=1"))

fit_allsites_flow_temp_mean <- readRDS(url("https://www.dropbox.com/scl/fi/uy6j279h4jvc5qbrl5iqb/fit_allsites_flow_temp_mean.rds?rlkey=a0qykz5ofu5of8zfqovoz3bca&dl=1"))
fit_allsites_flow_temp_cv <- readRDS(url("https://www.dropbox.com/scl/fi/l1zgikw4kbsf3sijsb9px/fit_allsites_flow_temp_cv.rds?rlkey=8kv6vs521vrqmsp48adn3sbpl&dl=1"))
fit_allsites_flow_water_mean <- readRDS(url("https://www.dropbox.com/scl/fi/2sxmqhnx926s76g5nizd6/fit_allsites_flow_water_mean.rds?rlkey=juhnekcasyd3c529epij5c0bl&dl=1"))
fit_allsites_flow_water_cv <- readRDS(url("https://www.dropbox.com/scl/fi/2nlxdyc3m8kkpbigjsnch/fit_allsites_flow_water_cv.rds?rlkey=5qsmx5bed8qxqtrikwjvr1evg&dl=1"))

fit_allsites_spi_temp_mean <- readRDS(url("https://www.dropbox.com/scl/fi/38ke36bs9hl4dlm7qcm1r/fit_allsites_spi_temp_mean.rds?rlkey=azsrlhwaj4bvh9cud720g24ip&dl=1"))
fit_allsites_spi_temp_cv <- readRDS(url("https://www.dropbox.com/scl/fi/rlgs4fedpxh18498sd743/fit_allsites_spi_temp_cv.rds?rlkey=wa5sygu748j07eybu34ud75og&dl=1"))
fit_allsites_spi_water_mean <- readRDS(url("https://www.dropbox.com/scl/fi/emsspmf96gab5281k6nw5/fit_allsites_spi_water_mean.rds?rlkey=qqwrevabchab1ao8clq3kqxes&dl=1"))
fit_allsites_spi_water_cv <- readRDS(url("https://www.dropbox.com/scl/fi/tpyjjtbyzqdq4xrrefzlb/fit_allsites_spi_water_cv.rds?rlkey=hsb4zlrhsm1m6vx7crlwhv0ly&dl=1"))

fit_allsites_grow_temp_mean <- readRDS(url("https://www.dropbox.com/scl/fi/s2og2sjw0ixvlgg2mib8p/fit_allsites_grow_temp_mean.rds?rlkey=diffkrgr3oj6ow01c8ky6aol2&dl=1"))
fit_allsites_grow_temp_cv <- readRDS(url("https://www.dropbox.com/scl/fi/ewoxlo2rhychdbpns56hm/fit_allsites_grow_temp_cv.rds?rlkey=b6e86w3k9uz2sdi29lsmx4brj&dl=1"))
fit_allsites_grow_water_mean <- readRDS(url("https://www.dropbox.com/scl/fi/afkdua2y5hzvvvmtxv1ch/fit_allsites_grow_water_mean.rds?rlkey=q49zprnlkkr3k0n9c4yhd4q9j&dl=1"))
fit_allsites_grow_water_cv <- readRDS(url("https://www.dropbox.com/scl/fi/qj19boa0mlr125cg55tlp/fit_allsites_grow_water_cv.rds?rlkey=8qynaksyib8gg91bjw9zaae5o&dl=1"))

## Chains convergence
mcmc_trace(fit_allsites_surv_temp_mean, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                                         bendotemp_s,bherbtemp_s,bendoherb_s))+theme_bw()
mcmc_pairs(fit_allsites_surv_temp_mean, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,bherbtemp_s,bendotemp_s,bendoherb_s))

mcmc_trace(fit_allsites_surv_temp_cv, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                                       bendotemp_s,bherbtemp_s,bendoherb_s))+theme_bw()
mcmc_pairs(fit_allsites_surv_temp_cv, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                        bendotemp_s,bherbtemp_s,bendoherb_s))
mcmc_trace(fit_allsites_surv_water_mean, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                                         bendotemp_s,bherbtemp_s,bendoherb_s))+theme_bw()
mcmc_pairs(fit_allsites_surv_water_mean, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                           bendotemp_s,bherbtemp_s,bendoherb_s))
mcmc_trace(fit_allsites_surv_water_cv, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                           bendotemp_s,bherbtemp_s,bendoherb_s))+theme_bw()
mcmc_pairs(fit_allsites_surv_water_cv, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                          bendotemp_s,bherbtemp_s,bendoherb_s))

mcmc_trace(fit_allsites_flow_temp_mean, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,btemp_f,
                                                                         bendotemp_f,bherbtemp_f,bendoherb_f))+theme_bw()
mcmc_trace(fit_allsites_flow_temp_cv, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,btemp_f,
                                                                       bendotemp_f,bherbtemp_f,bendoherb_f))+theme_bw()
mcmc_trace(fit_allsites_flow_water_mean, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,btemp_f,
                                                                          bendotemp_f,bherbtemp_f,bendoherb_f))+theme_bw()
mcmc_trace(fit_allsites_flow_water_cv, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,btemp_f,
                                                                        bendotemp_f,bherbtemp_f,bendoherb_f))+theme_bw()

mcmc_trace(fit_allsites_spi_temp_mean, pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,btemp_spk,
                                                                        bendotemp_spk,bherbtemp_spk,bendoherb_spk))+theme_bw()
mcmc_pairs(fit_allsites_spi_temp_mean, pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,btemp_spk,
                                                         bendotemp_spk,bherbtemp_spk,bendoherb_spk))
mcmc_trace(fit_allsites_spi_temp_cv, pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,btemp_spk,
                                                                      bendotemp_spk,bherbtemp_spk,bendoherb_spk))+theme_bw()
mcmc_pairs(fit_allsites_spi_temp_cv, pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,btemp_spk,
                                                       bendotemp_spk,bherbtemp_spk,bendoherb_spk))
mcmc_trace(fit_allsites_spi_water_mean,pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,btemp_spk,
                                                                         bendotemp_spk,bherbtemp_spk,bendoherb_spk))+theme_bw()
mcmc_pairs(fit_allsites_spi_water_mean,pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,btemp_spk,
                                                         bendotemp_spk,bherbtemp_spk,bendoherb_spk))
mcmc_trace(fit_allsites_spi_water_cv, pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,btemp_spk,
                                                                       bendotemp_spk,bherbtemp_spk,bendoherb_spk))+theme_bw()
mcmc_pairs(fit_allsites_spi_water_cv, pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,btemp_spk,
                                                        bendotemp_spk,bherbtemp_spk,bendoherb_spk))

mcmc_trace(fit_allsites_grow_temp_mean, pars = quote_bare(b0_g,bendo_g,bherb_g,btemp_g,btemp_g,bendotemp_g,
                                                                         bherbtemp_g,bendoherb_g))+theme_bw()
mcmc_trace(fit_allsites_grow_temp_cv, pars = quote_bare(b0_g,bendo_g,bherb_g,btemp_g,btemp_g,bendotemp_g,
                                                                       bherbtemp_g,bendoherb_g))+theme_bw()
mcmc_trace(fit_allsites_grow_water_mean, pars = quote_bare(b0_g,bendo_g,bherb_g,btemp_g,btemp_g,bendotemp_g,
                                                                          bherbtemp_g,bendoherb_g))+theme_bw()
mcmc_trace(fit_allsites_grow_water_cv, pars = quote_bare(b0_g,bendo_g,bherb_g,btemp_g,btemp_g,bendotemp_g,
                                                                        bherbtemp_g,bendoherb_g))+theme_bw()


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
phi_f_Tmean<- rstan::extract(fit_allsites_flow_temp_mean, pars = c("phi_f"))$phi_f
predF_Tcv <- rstan::extract(fit_allsites_flow_temp_cv, pars = c("predF"))$predF
phi_f_Tcv<- rstan::extract(fit_allsites_flow_temp_cv, pars = c("phi_f"))$phi_f
predF_Wmean <- rstan::extract(fit_allsites_flow_water_mean, pars = c("predF"))$predF
phi_f_Wmean<- rstan::extract(fit_allsites_flow_water_mean, pars = c("phi_f"))$phi_f
predF_Wcv <- rstan::extract(fit_allsites_flow_water_cv, pars = c("predF"))$predF
phi_f_Wcv<- rstan::extract(fit_allsites_flow_water_cv, pars = c("phi_f"))$phi_f

predSP_Tmean <- rstan::extract(fit_allsites_spi_temp_mean, pars = c("predF"))$predF
phi_spk_Tmean<- rstan::extract(fit_allsites_spi_temp_mean, pars = c("phi_spk"))$phi_spk
predSP_Tcv <- rstan::extract(fit_allsites_spi_temp_cv, pars = c("predF"))$predF
phi_spk_Tcv<- rstan::extract(fit_allsites_spi_temp_cv, pars = c("phi_spk"))$phi_spk
predSP_Wmean <- rstan::extract(fit_allsites_spi_water_mean, pars = c("predF"))$predF
phi_spk_Wmean<- rstan::extract(fit_allsites_spi_water_mean, pars = c("phi_spk"))$phi_spk
predSP_Wcv <- rstan::extract(fit_allsites_spi_water_cv, pars = c("predF"))$predF
phi_spk_Wmean<- rstan::extract(fit_allsites_spi_water_cv, pars = c("phi_spk"))$phi_spk

#draw 500 random samples from the joint posterior
n_post_draws <- 500
post_draws <- sample.int(dim(predS_Tmean)[1], n_post_draws)
#set up simulation output
y_s_tm_sim <- matrix(NA,n_post_draws,length(data_sites_surv_temp_mean$y_s))
y_s_tcv_sim <- matrix(NA,n_post_draws,length(data_sites_surv_temp_cv$y_s))
y_s_wm_sim <- matrix(NA,n_post_draws,length(data_sites_surv_water_mean$y_s))
y_s_wcv_sim <- matrix(NA,n_post_draws,length(data_sites_surv_water_cv$y_s))

y_g_tm_sim <- matrix(NA,n_post_draws,length(data_sites_grow_temp_mean$y_g))
y_g_tcv_sim <- matrix(NA,n_post_draws,length(data_sites_grow_temp_cv$y_g))
y_g_wm_sim <- matrix(NA,n_post_draws,length(data_sites_grow_water_mean$y_g))
y_g_wcv_sim <- matrix(NA,n_post_draws,length(data_sites_grow_water_cv$y_g))

y_f_tm_sim <- matrix(NA,n_post_draws,length(data_sites_flow_temp_mean$y_f))
y_f_tcv_sim <- matrix(NA,n_post_draws,length(data_sites_flow_temp_cv$y_f))
y_f_wm_sim <- matrix(NA,n_post_draws,length(data_sites_flow_water_mean$y_f))
y_f_wcv_sim <- matrix(NA,n_post_draws,length(data_sites_flow_water_cv$y_f))

y_sp_tm_sim <- matrix(NA,n_post_draws,length(data_sites_spi_temp_mean$y_spk))
y_sp_tcv_sim <- matrix(NA,n_post_draws,length(data_sites_spi_temp_cv$y_spk))
y_sp_wm_sim <- matrix(NA,n_post_draws,length(data_sites_spi_water_mean$y_spk))
y_sp_wcv_sim <- matrix(NA,n_post_draws,length(data_sites_spi_water_cv$y_spk))


#loop over the posterior and generate new observations
for(i in 1:n_post_draws){
  print(i)
  ## sample survival data (bernoulli)
  y_s_tm_sim[i,] <- rbinom(n=length(data_sites_surv_temp_mean$y_s), size=1, prob = invlogit(predS_Tmean[i,]))
  y_s_tcv_sim[i,] <- rbinom(n=length(data_sites_surv_temp_cv$y_s), size=1, prob = invlogit(predS_Tcv[i,]))
  y_s_wm_sim[i,] <- rbinom(n=length(data_sites_surv_water_mean$y_s), size=1, prob = invlogit(predS_Wmean[i,]))
  y_s_wcv_sim[i,] <- rbinom(n=length(data_sites_surv_water_cv$y_s), size=1, prob = invlogit(predS_Wcv[i,]))
  ## sample growth data (normal)
  y_g_tm_sim[i,] <- rnorm(n=length(data_sites_grow_temp_mean$y_g), mean = predG_Tmean[i,], sd=sigma_Tmean[i])
  y_g_tcv_sim[i,] <- rnorm(n=length(data_sites_grow_temp_cv$y_g), mean = predG_Tcv[i,], sd=sigma_Tcv[i])
  y_g_wm_sim[i,] <- rnorm(n=length(data_sites_grow_water_mean$y_g), mean = predG_Wmean[i,], sd=sigma_Wmean[i])
  y_g_wcv_sim[i,] <- rnorm(n=length(data_sites_grow_water_cv$y_g), mean = predG_Wcv[i,], sd=sigma_Wcv[i])
  ## sample flowering data (negative binomial)
  y_f_tm_sim[i,] <- rnbinom(n=length(data_sites_flow_temp_mean$y_f), mu = exp(predF_Tmean[i,]), size=phi_f_Tmean[i])
  y_f_tcv_sim[i,] <- rnbinom(n=length(data_sites_flow_temp_cv$y_f), mu = exp(predF_Tcv[i,]), size=phi_f_Tcv[i])
  y_f_wm_sim[i,] <- rnbinom(n=length(data_sites_flow_water_mean$y_f), mu = exp(predF_Wmean[i,]), size=phi_f_Wmean[i])
  y_f_wcv_sim[i,] <- rnbinom(n=length(data_sites_flow_water_cv$y_f), mu = exp(predF_Wcv[i,]), size=phi_f_Tcv[i])
  ## sample spikelet data (negative binomial)
  y_sp_tm_sim[i,] <- rnbinom(n=length(data_sites_spi_temp_mean$y_spk), mu = exp(predSP_Tmean[i,]), size=phi_spk_Tmean[i])
  y_sp_tcv_sim[i,] <- rnbinom(n=length(data_sites_spi_temp_cv$y_spk), mu = exp(predSP_Tcv[i,]), size=phi_spk_Tcv)
  y_sp_wm_sim[i,] <- rnbinom(n=length(data_sites_spi_water_mean$y_spk), mu = exp(predSP_Wmean[i,]), size=phi_spk_Wmean)
  y_sp_wcv_sim[i,] <- rnbinom(n=length(data_sites_spi_water_cv$y_spk), mu = exp(predSP_Wcv[i,]), size=phi_spk_Tcv)
}

#plot ppc overlay
#color_scheme_set("red")
color_scheme_set("blue")
bayesplot::ppc_dens_overlay(data_sites_surv_temp_mean$y_s, y_s_tm_sim)+
  xlab("Survival ")+
  ylab("Density")+
  ggtitle(("Mean temperature"))+
  theme_bw()+
  theme(legend.position = c(0.15, 0.7))->ppc_surv_temp_mean
ppc_dens_overlay(data_sites_surv_temp_mean$y_s, y_s_tm_sim)+
  xlab("Survival ")+
  ylab("Density")+
  ggtitle(("CV temperature"))+
  theme_bw()+
  theme(legend.position = "none")->ppc_surv_temp_cv
ppc_dens_overlay(data_sites_surv_temp_mean$y_s, y_s_tm_sim)+
  xlab("Survival ")+
  ylab("Density")+
  ggtitle(("Mean soil moisture"))+
  theme_bw()+
  theme(legend.position = "none")->ppc_surv_water_mean
ppc_dens_overlay(data_sites_surv_temp_mean$y_s, y_s_tm_sim)+
  xlab("Survival ")+
  ylab("Density")+
  ggtitle(("CV soil moisture"))+
  theme_bw()+
  theme(legend.position = "none")->ppc_surv_water_cv
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PPC_survival.pdf",useDingbats = F,height=5,width=6)
multiplot(ppc_surv_temp_mean,ppc_surv_temp_cv,ppc_surv_water_mean,
          ppc_surv_water_cv,cols=2)
dev.off()

bayesplot::ppc_dens_overlay(data_sites_grow_temp_mean$y_g, y_g_tm_sim)+
  xlab("Growth ")+
  ylab("Density")+
  ggtitle(("Mean temperature"))+
  theme_bw()+
  theme(legend.position = c(0.15, 0.7))->ppc_grow_temp_mean
ppc_dens_overlay(data_sites_grow_temp_mean$y_g, y_g_tm_sim)+
  xlab("Growth ")+
  ylab("Density")+
  ggtitle(("CV temperature"))+
  theme_bw()+
  theme(legend.position = "none")->ppc_grow_temp_cv
ppc_dens_overlay(data_sites_grow_temp_mean$y_g, y_g_tm_sim)+
  xlab("Growth ")+
  ylab("Density")+
  ggtitle(("Mean soil moisture"))+
  theme_bw()+
  theme(legend.position = "none")->ppc_grow_water_mean
ppc_dens_overlay(data_sites_grow_temp_mean$y_g, y_g_tm_sim)+
  xlab("Growth ")+
  ylab("Density")+
  ggtitle(("CV soil moisture"))+
  theme_bw()+
  theme(legend.position = "none")->ppc_grow_water_cv
# pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PPC_growth.pdf",useDingbats = F,height=5,width=6)
# multiplot(ppc_grow_temp_mean,ppc_grow_temp_cv,ppc_grow_water_mean,
#           ppc_grow_water_cv,cols=2)
# dev.off()

bayesplot::ppc_dens_overlay(data_sites_flow_temp_mean$y_f, y_f_tm_sim)+
  xlab("Flowering")+
  ylab("Density")+
  xlim(-10, 40)+
  ggtitle(("Mean temperature"))+
  theme_bw()+
  theme(legend.position = c(0.8, 0.7))->ppc_flow_temp_mean
ppc_dens_overlay(data_sites_flow_temp_mean$y_f, y_f_tm_sim)+
  xlab("Flowering")+
  ylab("Density")+
  xlim(-10, 40)+
  ggtitle(("CV temperature"))+
  theme_bw()+
  theme(legend.position = "none")->ppc_flow_temp_cv
ppc_dens_overlay(data_sites_flow_temp_mean$y_f, y_f_tm_sim)+
  xlab("Flowering")+
  ylab("Density")+
  xlim(-10, 40)+
  ggtitle(("Mean soil moisture"))+
  theme_bw()+
  theme(legend.position = "none")->ppc_flow_water_mean
ppc_dens_overlay(data_sites_flow_temp_mean$y_f, y_f_tm_sim)+
  xlab("Flowering")+
  ylab("Density")+
  xlim(-10, 40)+
  ggtitle(("CV soil moisture"))+
  theme_bw()+
  theme(legend.position = "none")->ppc_flow_water_cv
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PPC_flowering.pdf",useDingbats = F,height=5,width=6)
multiplot(ppc_flow_temp_mean,ppc_flow_temp_cv,ppc_flow_water_mean,
          ppc_flow_water_cv,cols=2)
dev.off()

bayesplot::ppc_dens_overlay(data_sites_spi_temp_mean$y_spk, y_sp_tm_sim)+
  xlab("Spikelet")+
  ylab("Density")+
  xlim(-10, 200)+
  ggtitle(("Mean temperature"))+
  theme_bw()+
  theme(legend.position = c(0.8, 0.7))->ppc_spi_temp_mean
ppc_dens_overlay(data_sites_spi_temp_mean$y_spk, y_sp_tm_sim)+
  xlab("Spikelet")+
  ylab("Density")+
  xlim(-10, 200)+
  ggtitle(("CV temperature"))+
  theme_bw()+
  theme(legend.position = "none")->ppc_spi_temp_cv
ppc_dens_overlay(data_sites_spi_temp_mean$y_spk, y_sp_tm_sim)+
  xlab("Spikelet")+
  ylab("Density")+
  xlim(-10, 200)+
  ggtitle(("Mean soil moisture"))+
  theme_bw()+
  theme(legend.position = "none")->ppc_spi_water_mean
ppc_dens_overlay(data_sites_spi_temp_mean$y_spk, y_sp_tm_sim)+
  xlab("Spikelet")+
  ylab("Density")+
  xlim(-10, 200)+
  ggtitle(("CV soil moisture"))+
  theme_bw()+
  theme(legend.position = "none")->ppc_spi_water_cv
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PPC_spikelet.pdf",useDingbats = F,height=5,width=6)
multiplot(ppc_spi_temp_mean,ppc_spi_temp_cv,ppc_spi_water_mean,
          ppc_spi_water_cv,cols=2)
dev.off()


## Model selection -----
## Model comparison based on epdl/looic 
### Survival
log_lik_surv_temp_mean <- loo::extract_log_lik(fit_allsites_surv_temp_mean, merge_chains = FALSE)
r_eff_surv_temp_mean <- loo::relative_eff(exp(log_lik_surv_temp_mean))
loo_surv_temp_mean <- loo(log_lik_surv_temp_mean, r_eff = r_eff_surv_temp_mean, cores = 4)
#plot(loo_surv_temp_mean)

log_lik_surv_temp_cv <- loo::extract_log_lik(fit_allsites_surv_temp_cv, merge_chains = FALSE)
r_eff_surv_temp_cv <- loo::relative_eff(exp(log_lik_surv_temp_cv))
loo_surv_temp_cv <- loo(log_lik_surv_temp_cv, r_eff = r_eff_surv_temp_cv, cores = 4)
#plot(loo_surv_temp_cv)

log_lik_surv_water_mean <- loo::extract_log_lik(fit_allsites_surv_water_mean, merge_chains = FALSE)
r_eff_surv_water_mean <- loo::relative_eff(exp(log_lik_surv_water_mean))
loo_surv_water_mean <- loo(log_lik_surv_water_mean, r_eff = r_eff_surv_water_mean, cores = 4)
#plot(loo_surv_water_mean)

log_lik_surv_water_cv <- loo::extract_log_lik(fit_allsites_surv_water_cv, merge_chains = FALSE)
r_eff_surv_water_cv <- loo::relative_eff(exp(log_lik_surv_water_cv))
loo_surv_water_cv <- loo(log_lik_surv_water_cv, r_eff = r_eff_surv_water_cv, cores = 4)
#plot(loo_surv_water_cv)

(comp_surv <- loo::loo_compare(loo_surv_temp_mean,loo_surv_temp_cv, loo_surv_water_mean,loo_surv_water_cv))

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

(comp_grow <- loo::loo_compare(loo_grow_temp_mean,loo_grow_temp_cv, loo_grow_water_mean,loo_grow_water_cv))

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

(comp_flow <- loo::loo_compare(loo_flow_temp_mean,loo_flow_temp_cv, loo_flow_water_mean,loo_flow_water_cv))

## Spiekelet
log_lik_spi_temp_mean <- loo::extract_log_lik(fit_allsites_spi_temp_mean, merge_chains = FALSE)
r_eff_spi_temp_mean <- loo::relative_eff(exp(log_lik_spi_temp_mean))
loo_spi_temp_mean <- loo(log_lik_spi_temp_mean, r_eff = r_eff_spi_temp_mean, cores = 4)
#plot(loo_spi_temp_mean)

log_lik_spi_temp_cv <- loo::extract_log_lik(fit_allsites_spi_temp_cv, merge_chains = FALSE)
r_eff_spi_temp_cv <- loo::relative_eff(exp(log_lik_spi_temp_cv))
loo_spi_temp_cv <- loo(log_lik_spi_temp_cv, r_eff = r_eff_spi_temp_cv, cores = 4)
#plot(loo_spi_temp_cv)

log_lik_spi_water_mean <- loo::extract_log_lik(fit_allsites_spi_water_mean, merge_chains = FALSE)
r_eff_spi_water_mean <- loo::relative_eff(exp(log_lik_spi_water_mean))
loo_spi_water_mean <- loo(log_lik_spi_water_mean, r_eff = r_eff_spi_water_mean, cores = 4)
#plot(loo_spi_water_mean)

log_lik_spi_water_cv <- loo::extract_log_lik(fit_allsites_spi_water_cv, merge_chains = FALSE)
r_eff_spi_water_cv <- loo::relative_eff(exp(log_lik_spi_water_cv))
loo_spi_water_cv <- loo(log_lik_spi_water_cv, r_eff = r_eff_spi_water_cv, cores = 4)
#plot(loo_spi_water_cv)

(comp_spi <- loo::loo_compare(loo_spi_temp_mean,loo_spi_temp_cv, loo_spi_water_mean,loo_spi_water_cv))

## # Posterior mean values for each vital rate----
## Temp CV
posterior_surv <- as.array(fit_allsites_surv_temp_cv)
#color_scheme_set("red")
surv<-mcmc_intervals(posterior_surv, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,
                                                       bendotemp_s,bherbtemp_s,bendoherb_s)) + 
  ggplot2::scale_y_discrete(limits = c("b0_s","bendo_s","bherb_s","btemp_s",
                                       "bendotemp_s","bherbtemp_s","bendoherb_s"),
                            labels=c("b0_s"="Grand Mean",
                                     "bendo_s"="Endophyte",
                                     "bherb_s"="Herbivory",
                                     "bsizesex_s"="size:sex",
                                     "btemp_s"="Temperature",
                                     "bendotemp_s"="Endo:Temp",
                                     "bherbtemp_s"="Herb:Temp",
                                     "bendoherb_s"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Survival)")+
  xlim(-8,8)+
  ggtitle("A") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

posterior_grow <- as.array(fit_allsites_grow_temp_cv)
grow<-mcmc_intervals(posterior_grow, pars = quote_bare(b0_g,bendo_g,bherb_g,btemp_g,
                                                       bendotemp_g,bherbtemp_g,bendoherb_g)) + 
  ggplot2::scale_y_discrete(limits = c("b0_g","bendo_g","bherb_g","btemp_g",
                                       "bendotemp_g","bherbtemp_g","bendoherb_g"),
                            labels=c("b0_g"="Grand Mean",
                                     "bendo_g"="Endophyte",
                                     "bherb_g"="Herbivory",
                                     "bsizesex_g"="size:sex",
                                     "btemp_g"="Temperature",
                                     "bendotemp_g"="Endo:Temp",
                                     "bherbtemp_g"="Herb:Temp",
                                     "bendoherb_g"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Growth)")+
  xlim(-2,2)+
  ggtitle("B") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

posterior_flow <- as.array(fit_allsites_flow_temp_cv)
flow<-mcmc_intervals(posterior_flow, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,
                                                       bendotemp_f,bherbtemp_f,bendoherb_f)) + 
  ggplot2::scale_y_discrete(limits = c("b0_f","bendo_f","bherb_f","btemp_f",
                                       "bendotemp_f","bherbtemp_f","bendoherb_f"),
                            labels=c("b0_f"="Grand Mean",
                                     "bendo_f"="Endophyte",
                                     "bherb_f"="Herbivory",
                                     "bsizesex_f"="size:sex",
                                     "btemp_f"="Temperature",
                                     "bendotemp_f"="Endo:Temp",
                                     "bherbtemp_f"="Herb:Temp",
                                     "bendoherb_f"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Inflorescence)")+
  xlim(-6,6)+
  ggtitle("C") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

posterior_flow <- as.array(fit_allsites_flow_temp_cv)

posterior_spi <- as.array(fit_allsites_spi_temp_cv)
spi<-mcmc_intervals(posterior_spi, pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,
                                                       bendotemp_spk,bherbtemp_spk,bendoherb_spk)) + 
  ggplot2::scale_y_discrete(limits = c("b0_spk","bendo_spk","bherb_spk","btemp_spk",
                                       "bendotemp_spk","bherbtemp_spk","bendoherb_spk"),
                            labels=c("b0_spk"="Grand Mean",
                                     "bendo_spk"="Endophyte",
                                     "bherb_spk"="Herbivory",
                                     "bsizesex_spk"="size:sex",
                                     "btemp_spk"="Temperature",
                                     "bendotemp_spk"="Endo:Temp",
                                     "bherbtemp_spk"="Herb:Temp",
                                     "bendoherb_spk"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Spikelet)")+
  xlim(-5,5)+
  ggtitle("D") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Posterior_mean.pdf",useDingbats = F,height=7,width=8)
ggarrange(surv,grow,flow,spi + rremove("ylab"), ncol = 2, nrow = 2)
dev.off()

## Temp mean
posterior_surv <- as.array(fit_allsites_surv_temp_mean)
#color_scheme_set("red")
surv<-mcmc_intervals(posterior_surv, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,
                                                       bendotemp_s,bherbtemp_s,bendoherb_s)) + 
  ggplot2::scale_y_discrete(limits = c("b0_s","bendo_s","bherb_s","btemp_s",
                                       "bendotemp_s","bherbtemp_s","bendoherb_s"),
                            labels=c("b0_s"="Grand Mean",
                                     "bendo_s"="Endophyte",
                                     "bherb_s"="Herbivory",
                                     "bsizesex_s"="size:sex",
                                     "btemp_s"="Temperature",
                                     "bendotemp_s"="Endo:Temp",
                                     "bherbtemp_s"="Herb:Temp",
                                     "bendoherb_s"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Survival)")+
  xlim(-8,8)+
  ggtitle("A") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

posterior_grow <- as.array(fit_allsites_grow_temp_mean)
grow<-mcmc_intervals(posterior_grow, pars = quote_bare(b0_g,bendo_g,bherb_g,btemp_g,
                                                       bendotemp_g,bherbtemp_g,bendoherb_g)) + 
  ggplot2::scale_y_discrete(limits = c("b0_g","bendo_g","bherb_g","btemp_g",
                                       "bendotemp_g","bherbtemp_g","bendoherb_g"),
                            labels=c("b0_g"="Grand Mean",
                                     "bendo_g"="Endophyte",
                                     "bherb_g"="Herbivory",
                                     "bsizesex_g"="size:sex",
                                     "btemp_g"="Temperature",
                                     "bendotemp_g"="Endo:Temp",
                                     "bherbtemp_g"="Herb:Temp",
                                     "bendoherb_g"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Growth)")+
  xlim(-2,2)+
  ggtitle("B") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

posterior_flow <- as.array(fit_allsites_flow_temp_mean)
flow<-mcmc_intervals(posterior_flow, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,
                                                       bendotemp_f,bherbtemp_f,bendoherb_f)) + 
  ggplot2::scale_y_discrete(limits = c("b0_f","bendo_f","bherb_f","btemp_f",
                                       "bendotemp_f","bherbtemp_f","bendoherb_f"),
                            labels=c("b0_f"="Grand Mean",
                                     "bendo_f"="Endophyte",
                                     "bherb_f"="Herbivory",
                                     "bsizesex_f"="size:sex",
                                     "btemp_f"="Temperature",
                                     "bendotemp_f"="Endo:Temp",
                                     "bherbtemp_f"="Herb:Temp",
                                     "bendoherb_f"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Inflorescence)")+
  xlim(-6,6)+
  ggtitle("C") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

posterior_spi <- as.array(fit_allsites_spi_temp_mean)
spi<-mcmc_intervals(posterior_spi, pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,
                                                     bendotemp_spk,bherbtemp_spk,bendoherb_spk)) + 
  ggplot2::scale_y_discrete(limits = c("b0_spk","bendo_spk","bherb_spk","btemp_spk",
                                       "bendotemp_spk","bherbtemp_spk","bendoherb_spk"),
                            labels=c("b0_spk"="Grand Mean",
                                     "bendo_spk"="Endophyte",
                                     "bherb_spk"="Herbivory",
                                     "bsizesex_spk"="size:sex",
                                     "btemp_spk"="Temperature",
                                     "bendotemp_spk"="Endo:Temp",
                                     "bherbtemp_spk"="Herb:Temp",
                                     "bendoherb_spk"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Spikelet)")+
  xlim(-5,5)+
  ggtitle("D") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Posterior_mean_temp_mean.pdf",useDingbats = F,height=7,width=8)
ggarrange(surv,grow,flow,spi + rremove("ylab"), ncol = 2, nrow = 2)
dev.off()

## Water CV
posterior_surv <- as.array(fit_allsites_surv_water_cv)
#color_scheme_set("red")
surv<-mcmc_intervals(posterior_surv, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,
                                                       bendotemp_s,bherbtemp_s,bendoherb_s)) + 
  ggplot2::scale_y_discrete(limits = c("b0_s","bendo_s","bherb_s","btemp_s",
                                       "bendotemp_s","bherbtemp_s","bendoherb_s"),
                            labels=c("b0_s"="Grand Mean",
                                     "bendo_s"="Endophyte",
                                     "bherb_s"="Herbivory",
                                     "bsizesex_s"="size:sex",
                                     "btemp_s"="Temperature",
                                     "bendotemp_s"="Endo:Temp",
                                     "bherbtemp_s"="Herb:Temp",
                                     "bendoherb_s"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Survival)")+
  xlim(-8,8)+
  ggtitle("A") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

posterior_grow <- as.array(fit_allsites_grow_water_cv)
grow<-mcmc_intervals(posterior_grow, pars = quote_bare(b0_g,bendo_g,bherb_g,btemp_g,
                                                       bendotemp_g,bherbtemp_g,bendoherb_g)) + 
  ggplot2::scale_y_discrete(limits = c("b0_g","bendo_g","bherb_g","btemp_g",
                                       "bendotemp_g","bherbtemp_g","bendoherb_g"),
                            labels=c("b0_g"="Grand Mean",
                                     "bendo_g"="Endophyte",
                                     "bherb_g"="Herbivory",
                                     "bsizesex_g"="size:sex",
                                     "btemp_g"="Temperature",
                                     "bendotemp_g"="Endo:Temp",
                                     "bherbtemp_g"="Herb:Temp",
                                     "bendoherb_g"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Growth)")+
  xlim(-2,2)+
  ggtitle("B") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

posterior_flow <- as.array(fit_allsites_flow_water_cv)
flow<-mcmc_intervals(posterior_flow, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,
                                                       bendotemp_f,bherbtemp_f,bendoherb_f)) + 
  ggplot2::scale_y_discrete(limits = c("b0_f","bendo_f","bherb_f","btemp_f",
                                       "bendotemp_f","bherbtemp_f","bendoherb_f"),
                            labels=c("b0_f"="Grand Mean",
                                     "bendo_f"="Endophyte",
                                     "bherb_f"="Herbivory",
                                     "bsizesex_f"="size:sex",
                                     "btemp_f"="Temperature",
                                     "bendotemp_f"="Endo:Temp",
                                     "bherbtemp_f"="Herb:Temp",
                                     "bendoherb_f"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Inflorescence)")+
  xlim(-6,6)+
  ggtitle("C") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

posterior_spi <- as.array(fit_allsites_spi_water_cv)
spi<-mcmc_intervals(posterior_spi, pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,
                                                     bendotemp_spk,bherbtemp_spk,bendoherb_spk)) + 
  ggplot2::scale_y_discrete(limits = c("b0_spk","bendo_spk","bherb_spk","btemp_spk",
                                       "bendotemp_spk","bherbtemp_spk","bendoherb_spk"),
                            labels=c("b0_spk"="Grand Mean",
                                     "bendo_spk"="Endophyte",
                                     "bherb_spk"="Herbivory",
                                     "bsizesex_spk"="size:sex",
                                     "btemp_spk"="Temperature",
                                     "bendotemp_spk"="Endo:Temp",
                                     "bherbtemp_spk"="Herb:Temp",
                                     "bendoherb_spk"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Spikelet)")+
  xlim(-5,5)+
  ggtitle("D") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Posterior_mean_water_cv.pdf",useDingbats = F,height=7,width=8)
ggarrange(surv,grow,flow,spi + rremove("ylab"), ncol = 2, nrow = 2)
dev.off()

## Water mean
posterior_surv <- as.array(fit_allsites_surv_water_mean)
#color_scheme_set("red")
surv<-mcmc_intervals(posterior_surv, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,
                                                       bendotemp_s,bherbtemp_s,bendoherb_s)) + 
  ggplot2::scale_y_discrete(limits = c("b0_s","bendo_s","bherb_s","btemp_s",
                                       "bendotemp_s","bherbtemp_s","bendoherb_s"),
                            labels=c("b0_s"="Grand Mean",
                                     "bendo_s"="Endophyte",
                                     "bherb_s"="Herbivory",
                                     "bsizesex_s"="size:sex",
                                     "btemp_s"="Temperature",
                                     "bendotemp_s"="Endo:Temp",
                                     "bherbtemp_s"="Herb:Temp",
                                     "bendoherb_s"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Survival)")+
  xlim(-8,8)+
  ggtitle("A") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

posterior_grow <- as.array(fit_allsites_grow_water_mean)
grow<-mcmc_intervals(posterior_grow, pars = quote_bare(b0_g,bendo_g,bherb_g,btemp_g,
                                                       bendotemp_g,bherbtemp_g,bendoherb_g)) + 
  ggplot2::scale_y_discrete(limits = c("b0_g","bendo_g","bherb_g","btemp_g",
                                       "bendotemp_g","bherbtemp_g","bendoherb_g"),
                            labels=c("b0_g"="Grand Mean",
                                     "bendo_g"="Endophyte",
                                     "bherb_g"="Herbivory",
                                     "bsizesex_g"="size:sex",
                                     "btemp_g"="Temperature",
                                     "bendotemp_g"="Endo:Temp",
                                     "bherbtemp_g"="Herb:Temp",
                                     "bendoherb_g"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Growth)")+
  xlim(-2,2)+
  ggtitle("B") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

posterior_flow <- as.array(fit_allsites_flow_water_mean)
flow<-mcmc_intervals(posterior_flow, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,
                                                       bendotemp_f,bherbtemp_f,bendoherb_f)) + 
  ggplot2::scale_y_discrete(limits = c("b0_f","bendo_f","bherb_f","btemp_f",
                                       "bendotemp_f","bherbtemp_f","bendoherb_f"),
                            labels=c("b0_f"="Grand Mean",
                                     "bendo_f"="Endophyte",
                                     "bherb_f"="Herbivory",
                                     "bsizesex_f"="size:sex",
                                     "btemp_f"="Temperature",
                                     "bendotemp_f"="Endo:Temp",
                                     "bherbtemp_f"="Herb:Temp",
                                     "bendoherb_f"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Inflorescence)")+
  xlim(-6,6)+
  ggtitle("C") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

posterior_spi <- as.array(fit_allsites_spi_water_mean)
spi<-mcmc_intervals(posterior_spi, pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,
                                                     bendotemp_spk,bherbtemp_spk,bendoherb_spk)) + 
  ggplot2::scale_y_discrete(limits = c("b0_spk","bendo_spk","bherb_spk","btemp_spk",
                                       "bendotemp_spk","bherbtemp_spk","bendoherb_spk"),
                            labels=c("b0_spk"="Grand Mean",
                                     "bendo_spk"="Endophyte",
                                     "bherb_spk"="Herbivory",
                                     "bsizesex_spk"="size:sex",
                                     "btemp_spk"="Temperature",
                                     "bendotemp_spk"="Endo:Temp",
                                     "bherbtemp_spk"="Herb:Temp",
                                     "bendoherb_spk"="Endo:Herb"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Spikelet)")+
  xlim(-5,5)+
  ggtitle("D") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Posterior_mean_water_mean.pdf",useDingbats = F,height=7,width=8)
ggarrange(surv,grow,flow,spi + rremove("ylab"), ncol = 2, nrow = 2)
dev.off()

## re-format data for plotting
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
# Temperature
## Survival
demography_climate_elvi_surv %>% 
  filter(Herbivory=="1")->elvi_surv_no_herb
elvi_surv_no_herb_plot <- data_summary(elvi_surv_no_herb, varname="surv_t1", 
                                       groupnames=c("Endo", "temp_mean"))
# Convert dose to a factor variable
elvi_surv_no_herb_plot$temp_mean=round(elvi_surv_no_herb_plot$temp_mean,2)
elvi_surv_no_herb_plot$temp_mean=as.factor(elvi_surv_no_herb_plot$temp_mean)
elvi_surv_no_herb_plot$Endo=as.factor(elvi_surv_no_herb_plot$Endo)

demography_climate_elvi_surv %>% 
  filter(Herbivory=="2")->elvi_surv_herb
elvi_surv_herb_plot <- data_summary(elvi_surv_herb, varname="surv_t1", 
                                       groupnames=c("Endo", "temp_mean"))
# Convert dose to a factor variable
elvi_surv_herb_plot$temp_mean=round(elvi_surv_herb_plot$temp_mean,2)
elvi_surv_herb_plot$temp_mean=as.factor(elvi_surv_herb_plot$temp_mean)
elvi_surv_herb_plot$Endo=as.factor(elvi_surv_herb_plot$Endo)

sur_h<- ggplot(elvi_surv_no_herb_plot, aes(x=temp_mean, y=surv_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=surv_t1-sd, ymax=surv_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Temperature (C)", y = "Survival")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = c(0.2, 0.2),
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))

surv<- ggplot(elvi_surv_herb_plot, aes(x=temp_mean, y=surv_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=surv_t1-sd, ymax=surv_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Temperature (C)", y = "Growth")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))


## Growth
demography_climate_elvi_grow %>% 
  filter(Herbivory=="1")->elvi_grow_no_herb

demography_climate_elvi_grow %>% 
  filter(Herbivory=="2")->elvi_grow_herb

elvi_grow_no_herb_plot <- data_summary(elvi_grow_no_herb, varname="grow", 
                    groupnames=c("Endo", "temp_mean"))
# Convert dose to a factor variable
elvi_grow_no_herb_plot$temp_mean=round(elvi_grow_no_herb_plot$temp_mean,2)
elvi_grow_no_herb_plot$temp_mean=as.factor(elvi_grow_no_herb_plot$temp_mean)
elvi_grow_no_herb_plot$Endo=as.factor(elvi_grow_no_herb_plot$Endo)
elvi_grow_herb_plot <- data_summary(elvi_grow_herb, varname="grow", 
                    groupnames=c("Endo", "temp_mean"))
# Convert dose to a factor variable
elvi_grow_herb_plot$temp_mean=round(elvi_grow_herb_plot$temp_mean,2)
elvi_grow_herb_plot$temp_mean=as.factor(elvi_grow_herb_plot$temp_mean)
elvi_grow_herb_plot$Endo=as.factor(elvi_grow_herb_plot$Endo)
grow_h<- ggplot(elvi_grow_no_herb_plot, aes(x=temp_mean, y=grow, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=grow-sd, ymax=grow+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="Fenced", x="Temperature (C)", y = "Growth")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = c(0.2, 0.23),
        plot.title = element_text(face="bold",hjust = 0.5),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "grey90"))

grow<- ggplot(elvi_grow_herb_plot, aes(x=temp_mean, y=grow, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=grow-sd, ymax=grow+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="Unfenced", x="Temperature (C)", y = "Growth")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        plot.title = element_text(face="bold",hjust = 0.5),
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))

## inflorescence
demography_climate_elvi_flowering %>% 
  filter(Herbivory=="1")->elvi_flow_no_herb
demography_climate_elvi_flowering %>% 
  filter(Herbivory=="2")->elvi_flow_herb

elvi_flow_no_herb_plot <- data_summary(elvi_flow_no_herb, varname="flow_t1", 
                                       groupnames=c("Endo", "temp_mean"))
# Convert dose to a factor variable
elvi_flow_no_herb_plot$temp_mean=round(elvi_flow_no_herb_plot$temp_mean,2)
elvi_flow_no_herb_plot$temp_mean=as.factor(elvi_flow_no_herb_plot$temp_mean)
elvi_flow_no_herb_plot$Endo=as.factor(elvi_flow_no_herb_plot$Endo)
elvi_flow_herb_plot <- data_summary(elvi_flow_herb, varname="flow_t1", 
                                    groupnames=c("Endo", "temp_mean"))
# Convert dose to a factor variable
elvi_flow_herb_plot$temp_mean=round(elvi_flow_herb_plot$temp_mean,2)
elvi_flow_herb_plot$temp_mean=as.factor(elvi_flow_herb_plot$temp_mean)
elvi_flow_herb_plot$Endo=as.factor(elvi_flow_herb_plot$Endo)

flow_h<- ggplot(elvi_flow_no_herb_plot, aes(x=temp_mean, y=flow_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=flow_t1-sd, ymax=flow_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Temperature (C)", y = "# Inflorescences")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray95"))

flow<- ggplot(elvi_flow_herb_plot, aes(x=temp_mean, y=flow_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=flow_t1-sd, ymax=flow_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Temperature (C)", y = " # Inflorescences")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))

## spikelet
demography_climate_elvi_spikelet %>% 
  filter(Herbivory=="1")->elvi_spi_no_herb
demography_climate_elvi_spikelet %>% 
  filter(Herbivory=="2")->elvi_spi_herb

elvi_spi_no_herb_plot <- data_summary(elvi_spi_no_herb, varname="spikelet_t1", 
                                       groupnames=c("Endo", "temp_mean"))
# Convert dose to a factor variable
elvi_spi_no_herb_plot$temp_mean=round(elvi_spi_no_herb_plot$temp_mean,2)
elvi_spi_no_herb_plot$temp_mean=as.factor(elvi_spi_no_herb_plot$temp_mean)
elvi_spi_no_herb_plot$Endo=as.factor(elvi_spi_no_herb_plot$Endo)
elvi_spi_herb_plot <- data_summary(elvi_spi_herb, varname="spikelet_t1", 
                                    groupnames=c("Endo", "temp_mean"))
# Convert dose to a factor variable
elvi_spi_herb_plot$temp_mean=round(elvi_spi_herb_plot$temp_mean,2)
elvi_spi_herb_plot$temp_mean=as.factor(elvi_spi_herb_plot$temp_mean)
elvi_spi_herb_plot$Endo=as.factor(elvi_spi_herb_plot$Endo)

spi_h<- ggplot(elvi_spi_no_herb_plot, aes(x=temp_mean, y=spikelet_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=spikelet_t1-sd, ymax=spikelet_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Temperature (C)", y = "# Spikelets")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))

spi<- ggplot(elvi_spi_herb_plot, aes(x=temp_mean, y=spikelet_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=spikelet_t1-sd, ymax=spikelet_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Temperature (C)", y = "# Spikelets")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Fig_temp_mean.pdf",useDingbats = F,height=8.5,width=6)
ggarrange(grow_h, grow + rremove("ylab"),flow_h,flow+ rremove("ylab"),spi_h,spi + rremove("ylab"), ncol = 2, nrow = 3)
dev.off()


# Temperature CGV
## Survival
elvi_surv_no_herb_plot_t_cv <- data_summary(elvi_surv_no_herb, varname="surv_t1", 
                                       groupnames=c("Endo", "temp_cv"))
# Convert dose to a factor variable
elvi_surv_no_herb_plot_t_cv$temp_cv=round(elvi_surv_no_herb_plot_t_cv$temp_cv,2)
elvi_surv_no_herb_plot_t_cv$temp_cv=as.factor(elvi_surv_no_herb_plot_t_cv$temp_cv)
elvi_surv_no_herb_plot_t_cv$Endo=as.factor(elvi_surv_no_herb_plot_t_cv$Endo)

elvi_surv_herb_plot_t_cv <- data_summary(elvi_surv_herb, varname="surv_t1", 
                                    groupnames=c("Endo", "temp_cv"))
# Convert dose to a factor variable
elvi_surv_herb_plot_t_cv$temp_cv=round(elvi_surv_herb_plot_t_cv$temp_cv,2)
elvi_surv_herb_plot_t_cv$temp_cv=as.factor(elvi_surv_herb_plot_t_cv$temp_cv)
elvi_surv_herb_plot_t_cv$Endo=as.factor(elvi_surv_herb_plot_t_cv$Endo)

sur_h_t_cv<- ggplot(elvi_surv_no_herb_plot_t_cv, aes(x=temp_cv, y=surv_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=surv_t1-sd, ymax=surv_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Temperature (CV)", y = "Survival")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = c(0.2, 0.2),
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))

surv_t_cv<- ggplot(elvi_surv_herb_plot_t_cv, aes(x=temp_cv, y=surv_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=surv_t1-sd, ymax=surv_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Temperature (CV)", y = "Growth")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))


## Growth
elvi_grow_no_herb_plot_t_cv <- data_summary(elvi_grow_no_herb, varname="grow", 
                                       groupnames=c("Endo", "temp_cv"))
# Convert dose to a factor variable
elvi_grow_no_herb_plot_t_cv$temp_cv=round(elvi_grow_no_herb_plot_t_cv$temp_cv,2)
elvi_grow_no_herb_plot_t_cv$temp_cv=as.factor(elvi_grow_no_herb_plot_t_cv$temp_cv)
elvi_grow_no_herb_plot_t_cv$Endo=as.factor(elvi_grow_no_herb_plot_t_cv$Endo)

elvi_grow_herb_plot_t_cv <- data_summary(elvi_grow_herb, varname="grow", 
                                    groupnames=c("Endo", "temp_cv"))
# Convert dose to a factor variable
elvi_grow_herb_plot_t_cv$temp_cv=round(elvi_grow_herb_plot_t_cv$temp_cv,2)
elvi_grow_herb_plot_t_cv$temp_cv=as.factor(elvi_grow_herb_plot_t_cv$temp_cv)
elvi_grow_herb_plot_t_cv$Endo=as.factor(elvi_grow_herb_plot_t_cv$Endo)

grow_h_t_cv<- ggplot(elvi_grow_no_herb_plot_t_cv, aes(x=temp_cv, y=grow, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=grow-sd, ymax=grow+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="Fenced", x="Temperature (CV)", y = "Growth")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = c(0.2, 0.23),
        plot.title = element_text(face="bold",hjust = 0.5),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "grey90"))

grow_t_cv<- ggplot(elvi_grow_herb_plot_t_cv, aes(x=temp_cv, y=grow, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=grow-sd, ymax=grow+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="Unfenced", x="Temperature (CV)", y = "Growth")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        plot.title = element_text(face="bold",hjust = 0.5),
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))

## inflorescence
elvi_flow_no_herb_plot_t_cv <- data_summary(elvi_flow_no_herb, varname="flow_t1", 
                                       groupnames=c("Endo", "temp_cv"))
# Convert dose to a factor variable
elvi_flow_no_herb_plot_t_cv$temp_cv=round(elvi_flow_no_herb_plot_t_cv$temp_cv,2)
elvi_flow_no_herb_plot_t_cv$temp_cv=as.factor(elvi_flow_no_herb_plot_t_cv$temp_cv)
elvi_flow_no_herb_plot_t_cv$Endo=as.factor(elvi_flow_no_herb_plot_t_cv$Endo)

elvi_flow_herb_plot_t_cv <- data_summary(elvi_flow_herb, varname="flow_t1", 
                                    groupnames=c("Endo", "temp_cv"))
# Convert dose to a factor variable
elvi_flow_herb_plot_t_cv$temp_cv=round(elvi_flow_herb_plot_t_cv$temp_cv,2)
elvi_flow_herb_plot_t_cv$temp_cv=as.factor(elvi_flow_herb_plot_t_cv$temp_cv)
elvi_flow_herb_plot_t_cv$Endo=as.factor(elvi_flow_herb_plot_t_cv$Endo)

flow_h_cv<- ggplot(elvi_flow_herb_plot_t_cv, aes(x=temp_cv, y=flow_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=flow_t1-sd, ymax=flow_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Temperature (CV)", y = "# Inflorescences")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray95"))

flow_cv<- ggplot(elvi_flow_herb_plot_t_cv, aes(x=temp_cv, y=flow_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=flow_t1-sd, ymax=flow_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Temperature (CV)", y = " # Inflorescences")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))

## spikelet

elvi_spi_no_herb_plot_t_cv <- data_summary(elvi_spi_no_herb, varname="spikelet_t1", 
                                      groupnames=c("Endo", "temp_cv"))
# Convert dose to a factor variable
elvi_spi_no_herb_plot_t_cv$temp_cv=round(elvi_spi_no_herb_plot_t_cv$temp_cv,2)
elvi_spi_no_herb_plot_t_cv$temp_cv=as.factor(elvi_spi_no_herb_plot_t_cv$temp_cv)
elvi_spi_no_herb_plot_t_cv$Endo=as.factor(elvi_spi_no_herb_plot_t_cv$Endo)

elvi_spi_herb_plot_t_cv <- data_summary(elvi_spi_herb, varname="spikelet_t1", 
                                   groupnames=c("Endo", "temp_cv"))
# Convert dose to a factor variable
elvi_spi_herb_plot_t_cv$temp_cv=round(elvi_spi_herb_plot_t_cv$temp_cv,2)
elvi_spi_herb_plot_t_cv$temp_cv=as.factor(elvi_spi_herb_plot_t_cv$temp_cv)
elvi_spi_herb_plot_t_cv$Endo=as.factor(elvi_spi_herb_plot_t_cv$Endo)

spi_h_t_cv<- ggplot(elvi_spi_no_herb_plot_t_cv, aes(x=temp_cv, y=spikelet_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=spikelet_t1-sd, ymax=spikelet_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Temperature (CV)", y = "# Spikelets")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))

spi_t_cv<- ggplot(elvi_spi_herb_plot_t_cv, aes(x=temp_cv, y=spikelet_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=spikelet_t1-sd, ymax=spikelet_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Temperature (CV)", y = "# Spikelets")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Fig_temp_cv.pdf",useDingbats = F,height=8.5,width=6)
ggarrange(grow_h_t_cv, grow_t_cv + rremove("ylab"),flow_h_cv,flow_cv+ rremove("ylab"),spi_h_t_cv,spi_h_t_cv + rremove("ylab"), ncol = 2, nrow = 3)
dev.off()

# Soil moisture
## Survival
elvi_surv_no_herb_plot_w <- data_summary(elvi_surv_no_herb, varname="surv_t1", 
                                       groupnames=c("Endo", "water_mean"))
# Convert dose to a factor variable
elvi_surv_no_herb_plot_w$water_mean=round(elvi_surv_no_herb_plot_w$water_mean,2)
elvi_surv_no_herb_plot_w$water_mean=as.factor(elvi_surv_no_herb_plot_w$water_mean)
elvi_surv_no_herb_plot_w$Endo=as.factor(elvi_surv_no_herb_plot_w$Endo)

elvi_surv_herb_plot_w <- data_summary(elvi_surv_herb, varname="surv_t1", 
                                    groupnames=c("Endo", "water_mean"))
# Convert dose to a factor variable
elvi_surv_herb_plot_w$water_mean=round(elvi_surv_herb_plot_w$water_mean,2)
elvi_surv_herb_plot_w$water_mean=as.factor(elvi_surv_herb_plot_w$water_mean)
elvi_surv_herb_plot_w$Endo=as.factor(elvi_surv_herb_plot_w$Endo)

sur_h_w<- ggplot(elvi_surv_no_herb_plot_w, aes(x=water_mean, y=surv_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=surv_t1-sd, ymax=surv_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Soil moisture (%) ", y = "Survival")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = c(0.2, 0.2),
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))

surv_w<- ggplot(elvi_surv_herb_plot_w, aes(x=water_mean, y=surv_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=surv_t1-sd, ymax=surv_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Soil moisture (%)", y = "Growth")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))


## Growth
elvi_grow_no_herb_plot_w <- data_summary(elvi_grow_no_herb, varname="grow", 
                                       groupnames=c("Endo", "water_mean"))
# Convert dose to a factor variable
elvi_grow_no_herb_plot_w$water_mean=round(elvi_grow_no_herb_plot_w$water_mean,2)
elvi_grow_no_herb_plot_w$water_mean=as.factor(elvi_grow_no_herb_plot_w$water_mean)
elvi_grow_no_herb_plot_w$Endo=as.factor(elvi_grow_no_herb_plot_w$Endo)

elvi_grow_herb_plot_w <- data_summary(elvi_grow_herb, varname="grow", 
                                    groupnames=c("Endo", "water_mean"))
# Convert dose to a factor variable
elvi_grow_herb_plot_w$water_mean=round(elvi_grow_herb_plot_w$water_mean,2)
elvi_grow_herb_plot_w$water_mean=as.factor(elvi_grow_herb_plot_w$water_mean)
elvi_grow_herb_plot_w$Endo=as.factor(elvi_grow_herb_plot_w$Endo)
grow_h_w<- ggplot(elvi_grow_no_herb_plot_w, aes(x=water_mean, y=grow, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=grow-sd, ymax=grow+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="Fenced", x="TSoil moisture (%)", y = "Growth")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = c(0.4, 0.23),
        plot.title = element_text(face="bold",hjust = 0.5),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "grey90"))

grow_w<- ggplot(elvi_grow_herb_plot_w, aes(x=water_mean, y=grow, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=grow-sd, ymax=grow+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="Unfenced", x="Soil moisture (%)", y = "Growth")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        plot.title = element_text(face="bold",hjust = 0.5),
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))

## inflorescence
elvi_flow_no_herb_plot_w <- data_summary(elvi_flow_no_herb, varname="flow_t1", 
                                       groupnames=c("Endo", "water_mean"))
# Convert dose to a factor variable
elvi_flow_no_herb_plot_w$water_mean=round(elvi_flow_no_herb_plot_w$water_mean,2)
elvi_flow_no_herb_plot_w$water_mean=as.factor(elvi_flow_no_herb_plot_w$water_mean)
elvi_flow_no_herb_plot_w$Endo=as.factor(elvi_flow_no_herb_plot_w$Endo)

elvi_flow_herb_plot_w <- data_summary(elvi_flow_herb, varname="flow_t1", 
                                    groupnames=c("Endo", "water_mean"))
# Convert dose to a factor variable
elvi_flow_herb_plot_w$water_mean=round(elvi_flow_herb_plot_w$water_mean,2)
elvi_flow_herb_plot_w$water_mean=as.factor(elvi_flow_herb_plot_w$water_mean)
elvi_flow_herb_plot_w$Endo=as.factor(elvi_flow_herb_plot_w$Endo)

flow_h_w<- ggplot(elvi_flow_no_herb_plot_w, aes(x=water_mean, y=flow_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=flow_t1-sd, ymax=flow_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Soil moisture (%)", y = "# Inflorescences")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray95"))

flow_w<- ggplot(elvi_flow_herb_plot_w, aes(x=water_mean, y=flow_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=flow_t1-sd, ymax=flow_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Soil moisture (%)", y = " # Inflorescences")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))

## spikelet
elvi_spi_no_herb_plot_w <- data_summary(elvi_spi_no_herb, varname="spikelet_t1", 
                                      groupnames=c("Endo", "water_mean"))
# Convert dose to a factor variable
elvi_spi_no_herb_plot_w$water_mean=round(elvi_spi_no_herb_plot_w$water_mean,2)
elvi_spi_no_herb_plot_w$water_mean=as.factor(elvi_spi_no_herb_plot_w$water_mean)
elvi_spi_no_herb_plot_w$Endo=as.factor(elvi_spi_no_herb_plot_w$Endo)
elvi_spi_herb_plot_w <- data_summary(elvi_spi_herb, varname="spikelet_t1", 
                                   groupnames=c("Endo", "water_mean"))
# Convert dose to a factor variable
elvi_spi_herb_plot_w$water_mean=round(elvi_spi_herb_plot_w$water_mean,2)
elvi_spi_herb_plot_w$water_mean=as.factor(elvi_spi_herb_plot_w$water_mean)
elvi_spi_herb_plot_w$Endo=as.factor(elvi_spi_herb_plot_w$Endo)

spi_h_w<- ggplot(elvi_spi_no_herb_plot_w, aes(x=water_mean, y=spikelet_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=spikelet_t1-sd, ymax=spikelet_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Soil moisture (%)", y = "# Spikelets")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))

spi_w<- ggplot(elvi_spi_herb_plot_w, aes(x=water_mean, y=spikelet_t1, group=Endo, color=Endo)) + 
  geom_line() +
  geom_point(size = 3)+
  guides(color = guide_legend(title = "Endophyte")) +
  geom_errorbar(aes(ymin=spikelet_t1-sd, ymax=spikelet_t1+sd), width=.3,
                position=position_dodge(0.03))+
  labs(title="", x="Soil moisture (%)", y = "# Spikelets")+
  theme_bw() +
  scale_color_manual(name="Endophyte",labels=c("1" = expression(E^"-"), "2"=expression(E^"+")), values=c('#999999','#E69F00'))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", linetype = "blank",fill = "gray90"))

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Fig_water_mean.pdf",useDingbats = F,height=8.5,width=6)
ggarrange(grow_h_w, grow_w + rremove("ylab"),flow_h_w,flow_w + rremove("ylab"),spi_h_w,spi_w + rremove("ylab"), ncol = 2, nrow = 3)
dev.off()

