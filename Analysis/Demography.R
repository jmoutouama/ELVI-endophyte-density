# Project: 
# Purpose: Fit vital rate models to test the effect of grass-endophyte symbiosis and endophyte hyphal density on  vital rate models (survival, growth, flowering,fertility).
# Note: Raster files are too large to provide in public repository. They are stored on a local machine
# Authors: Jacob Moutouama
# Date last modified (Y-M-D): 2024-08-03

rm(list = ls())

# load packages
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
library(rgdal)
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

jacob_path<-"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density"
# tom_path<-"C:/Users/tm9/Dropbox/github/ELVI-endophyte-density" 
choose_path<-jacob_path
# Demographic data -----
# Merge the demographic census
datini<-read_csv(paste0(choose_path,"/Data/Initialdata.csv"))
dat23<-read_csv(paste0(choose_path,"/Data/census2023.csv"))
dat24<-read_csv(paste0(choose_path,"/Data/census2024.csv"))
datherbivory<-read_csv(paste0(choose_path,"/Data/herbivory.csv"))
unique(datini$Site)
unique(datini$dat23)
unique(datini$dat24)
# calculate the total spikelet for each census
dat23 %>% 
  mutate(spikelet_23=rowSums(across(Spikelet_A:Spikelet_C)))->dat23_spike
dat24 %>% 
  mutate(spikelet_24=rowSums(across(Spikelet_A:Spikelet_C)),Inf_24=rowSums(across(attachedInf_24:brokenInf_24)))->dat24_spike

## Merge the initial data with the 23 data and the 23 data with the 24 -----

dat2324 <- left_join(x = datini,y =dat24_spike,by=c("Tag_ID"))
names(dat2324)
unique(dat2324$Site)
dat2324 %>% 
  mutate(tiller_t=ini_Tiller,
         tiller_t1=Tiller_24,
         inf_t1=Inf_24,
         spikelet_t1=spikelet_24,
         tiller_Herb=tiller_herb_24) %>% 
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
                inf_t1,
                spikelet_t1,
                stroma,
                tiller_Herb)->dat2324_all_sp

#names(dat2324_all_sp)
## Merge the demographic data with the herbivory data -----
dat2324_all_sp_herb<-left_join(x=dat2324_all_sp,y=datherbivory,by=c("Site","Plot","Species"))# Merge the demographic data with the herbivory data
view(dat2324_all_sp_herb)
# head(dat2324_all_sp_herb)

## Find the starting and ending dates are correct
dat2324_all_sp_herb %>% 
  dplyr::select(Site,Species,date_23,date_24) %>% 
  group_by(Site,Species) %>% 
  unique()->dat2324_all_sp_herb_dates

view(dat2324_all_sp_herb_dates)

# HOBO data ----
## format date and separate year-month-day
list.files(path = paste0(choose_path,"/Data/HOBO data/"),  
           pattern = "*.xlsx", full.names = TRUE) %>% # Identify all excel files
  lapply(read_excel) %>%                              # Store all files in list
  bind_rows ->hobo_data_raw # get HOBO data

tidyr::separate(hobo_data_raw, "date",
                into = c('longdate', 'time'),
                sep= ' ') %>%
  tidyr::separate('longdate', # Separate the ‘longdate’ column into separate columns for month, day and year using the separate() function.
                  into = c('year','month', 'day'),
                  sep= '-',
                  remove = FALSE)->hobo_data_full 

## double check if the starting and ending dates are correct
hobo_data_full %>% 
  group_by(site) %>% 
  summarise(start=range(longdate)[1],
            end=range(longdate)[2],
            duration=as.Date(end)-as.Date(start))->hobo_dates 

## average over days to look at overall trend across sites
hobo_data_full %>% 
  group_by(longdate,site,day) %>% 
  summarise(daily_mean_moist=mean(water),daily_mean_temp=mean(temperature))->HOBO_daily

## Plot the daily trend for temperature and soil moisture from start to end
hobo_means<-HOBO_daily %>% 
  group_by(site) %>% 
  summarise(mean_temp=mean(daily_mean_temp),
            mean_moisture=mean(daily_mean_moist))

data_plotclim<-data.frame(site=c(HOBO_daily$site,HOBO_daily$site),daily_mean_clim=c(HOBO_daily$daily_mean_temp,HOBO_daily$daily_mean_moist),date=c(HOBO_daily$longdate,HOBO_daily$longdate),clim=c(rep("temp",nrow(HOBO_daily)),rep("water",nrow(HOBO_daily))))
  
site_names <- c("LAF"="Lafayette",
                 "HUN"="Huntville",
                 "BAS"="Bastrop",
                "COL"="College Station",
                "KER" ="Kerville",
                "BFL" ="Brackenridge",
                "SON"="Sonora")


  
figtempsite<-ggplot(HOBO_daily, aes(x=as.Date(longdate, format= "%Y - %m - %d"), y=daily_mean_temp))+
  geom_line(aes(colour=site))+
  ggtitle("a")+
  scale_fill_jco()+
  theme_bw()+ 
  theme(legend.position = "none",
        axis.text.x = element_text(size=4.5,color="black", angle=0))+
  labs( y="Daily temperature  (°C)", x="")+
  # facet_grid(~factor(site,levels=c("LAF", "HUN", "BAS", "COL" ,"KER" ,"BLF", "SON")))+
  facet_grid(~site,labeller = labeller(site=site_names))+
  geom_hline(data=hobo_means,aes(yintercept = mean_temp,colour=site))

figmoistsite<-ggplot(HOBO_daily, aes(x=as.Date(longdate, format= "%Y - %m - %d"), y=daily_mean_moist))+
  geom_line(aes(colour=site))+
  ggtitle("b")+
  scale_fill_jco()+
  theme_bw()+ 
  theme(legend.position = "none",
        axis.text.x = element_text(size=4.55, color="black",angle=0))+
  labs( y="Daily soil moisture (wfv)", x="Month")+
  facet_grid(~site,labeller = labeller(site=site_names))+
  geom_hline(data=hobo_means,aes(yintercept = mean_moisture,colour=site))

# pdf("paste0(choose_path/Figure/climatesite.pdf",height =5,width =12,useDingbats = F)
# (Figclimatesite<-ggpubr::ggarrange(figtempsite,figmoistsite,common.legend = FALSE,ncol = 1, nrow = 2))
# dev.off()

## Average over days

HOBO_daily %>%
  group_by(site)%>%
  summarise(water_mean = mean(daily_mean_moist),
            water_cv=sd(daily_mean_moist)/water_mean,
            temp_mean = mean(daily_mean_temp),
            temp_cv=sd(daily_mean_temp)/temp_mean) %>% 
  left_join(.,hobo_dates,by=c("site"))->HOBO_summary


HOBO_summary %>% 
  rename(Site=site)->HOBO_summary_clean
## Merge the demographic data with the climatic data -----
demography_climate<-left_join(x=demography,y=HOBO_summary_clean,by=c("Site"))# Merge the demographic data with the temperature data
# Subset only ELVI data -----
demography_climate %>% 
  filter(Species=="ELVI")->demography_climate_elvi
demography_climate_elvi$surv1<-1*(!is.na(demography_climate_elvi$tiller_t) & !is.na(demography_climate_elvi$tiller_t1))
demography_climate_elvi$site_plot<-interaction(demography_climate_elvi$Site,demography_climate_elvi$Plot)
demography_climate_elvi$grow<-(log(demography_climate_elvi$tiller_t1+1) - log(demography_climate_elvi$tiller_t+1))# Relative growth rate
# names(demography_climate_elvi)
# view(demography_climate_elvi)
# summary(demography_climate_elvi)

# H1: We hypothesized that stress associated with aridity and low precipitation would strengthen the plant-fungal mutualism, such that the fitness benefits of endophyte symbiosis are maximized at the range edge. ----

## Survival----
## Read and format survival data to build the model
demography_climate_elvi %>% 
  subset( tiller_t > 0 )%>%
  dplyr::select(year, Population, Site, Plot,site_plot, Endo, Herbivory,
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

sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 3, 
  chains = 3
)

fit_allsites_surv_temp_mean <- stan(
 file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_survival.stan",
 data = data_sites_surv_temp_mean,
 warmup = sim_pars$warmup,
 iter = sim_pars$iter,
 thin = sim_pars$thin,
 chains = sim_pars$chains)

fit_allsites_surv_temp_cv <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_survival.stan",
  data = data_sites_surv_temp_cv,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

fit_allsites_surv_water_mean <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_survival.stan",
  data = data_sites_surv_water_mean,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

fit_allsites_surv_water_cv <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_survival.stan",
  data = data_sites_surv_water_cv,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

## Save RDS file for further use
# saveRDS(fit_allsites_surv_temp_mean, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Survival/fit_allsites_surv_temp_mean.rds')
# saveRDS(fit_allsites_surv_temp_cv, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Survival/fit_allsites_surv_temp_cv.rds')
# saveRDS(fit_allsites_surv_water_mean, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Survival/fit_allsites_surv_water_mean.rds')
# saveRDS(fit_allsites_surv_water_cv, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Survival/fit_allsites_surv_water_cv.rds')

## Flowering----
demography_climate_elvi %>% 
  subset( tiller_t1 > 0 )%>%
  dplyr::select(year, Population, Site, Plot,site_plot, Endo, Herbivory,
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

fit_allsites_flow_temp_mean <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_flowering.stan",
  data = data_sites_flow_temp_mean,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

fit_allsites_flow_temp_cv <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_flowering.stan",
  data = data_sites_flow_temp_cv,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

fit_allsites_flow_water_mean <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_flowering.stan",
  data = data_sites_flow_water_mean,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

fit_allsites_flow_water_cv <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_flowering.stan",
  data = data_sites_flow_water_cv,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

# saveRDS(fit_allsites_flow_temp_mean, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Flowering/fit_allsites_flow_temp_mean.rds')
# saveRDS(fit_allsites_flow_temp_cv, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Flowering/fit_allsites_flow_temp_cv.rds')
# saveRDS(fit_allsites_flow_water_mean, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Flowering/fit_allsites_flow_water_mean.rds')
# saveRDS(fit_allsites_flow_water_cv, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Flowering/fit_allsites_flow_water_cv.rds')

## Spikelet----
demography_climate_elvi %>% 
  subset( spikelet_t1 > 0 & tiller_t1 > 0 )%>%
  dplyr::select(year, Population, Site, Plot,site_plot, Endo, Herbivory,
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

data_sites_spklow_water_cv <- list( n_sites    = demography_climate_elvi_spikelet$Site %>% n_distinct,
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

fit_allsites_spi_temp_mean <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_spilkelet.stan",
  data = data_sites_spi_temp_mean,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

fit_allsites_spi_temp_cv <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_spilkelet.stan",
  data = data_sites_spi_temp_cv,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

fit_allsites_spi_water_mean <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_spilkelet.stan",
  data = data_sites_spi_water_mean,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

fit_allsites_spi_water_cv <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_spilkelet.stan",
  data = data_sites_spi_water_mean,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

#Save RDS file for further use
# saveRDS(fit_allsites_spi_temp_mean, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Spikelet/fit_allsites_spi_temp_mean.rds')
# saveRDS(fit_allsites_spi_temp_cv, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Spikelet/fit_allsites_spi_temp_cv.rds')
# saveRDS(fit_allsites_spi_water_mean, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Spikelet/fit_allsites_spi_water_mean.rds')
# saveRDS(fit_allsites_spi_water_cv, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Spikelet/fit_allsites_spi_water_cv.rds')

## Growth----
## Read and format the growth data to build the model
#hist(demography_climate_elvi$grow)
demography_climate_elvi %>% 
  subset( tiller_t > 0 & tiller_t1 > 0)%>%
  dplyr::select(year, Population, Site, Plot,site_plot, Endo, Herbivory,
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
data_sites_surv_water_mean <- list( n_sites    = demography_climate_elvi_grow$Site %>% n_distinct,
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
data_sites_surv_water_cv <- list( n_sites    = demography_climate_elvi_grow$Site %>% n_distinct,
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


fit_allsites_grow_temp_mean <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_growth.stan",
  data = data_sites_grow_temp_mean,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

fit_allsites_grow_temp_cv <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_growth.stan",
  data = data_sites_grow_temp_cv,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

fit_allsites_grow_water_mean <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_growth.stan",
  data = data_sites_surv_water_mean,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

fit_allsites_grow_water_cv <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_growth.stan",
  data = data_sites_surv_water_cv,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

# Save RDS file for further use
# saveRDS(fit_allsites_grow_temp_mean, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Growth/fit_allsites_grow_temp_mean.rds')
# saveRDS(fit_allsites_grow_temp_cv, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Growth/fit_allsites_grow_temp_cv.rds')
# saveRDS(fit_allsites_grow_water_mean, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Growth/fit_allsites_grow_water_mean.rds')
# saveRDS(fit_allsites_grow_water_cv, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/Stan output/Growth/fit_allsites_grow_water_cv.rds')

## load stan output 
fit_allsites_surv_temp_mean <- readRDS(url("https://www.dropbox.com/scl/fi/g3spz9lff9xy3153o96vo/fit_allsites_surv_temp_mean.rds?rlkey=jnv64cw87w1wi9te8brvfh7df&dl=1"))
fit_allsites_surv_temp_cv <- readRDS(url("https://www.dropbox.com/scl/fi/vq8ryu7vgtp3yqs68ewff/fit_allsites_surv_temp_cv.rds?rlkey=i227tu61cp224hkjpgekox11r&dl=1"))
fit_allsites_surv_water_mean <- readRDS(url("https://www.dropbox.com/scl/fi/at18rlo9e53pqc1tmmund/fit_allsites_surv_water_mean.rds?rlkey=hvvbxzn9g1x4brlidxcm147n3&dl=1"))
fit_allsites_surv_water_cv <- readRDS(url("https://www.dropbox.com/scl/fi/wdnd6e2oc4zunmmhsvayn/fit_allsites_surv_water_cv.rds?rlkey=sbim89fp48gnsoagpk1r57t4k&dl=1"))

fit_allsites_flow_temp_mean <- readRDS(url("https://www.dropbox.com/scl/fi/rlmw1uc9nsbic0reiuju4/fit_allsites_flow_temp_mean.rds?rlkey=rs4289erkiphqik23q70mx83n&dl=1"))
fit_allsites_flow_temp_cv <- readRDS(url("https://www.dropbox.com/scl/fi/ba574pxx58fh03feo6c3m/fit_allsites_flow_temp_cv.rds?rlkey=dsnjc4b6pppe9ttp1sqbzqggn&dl=1"))
fit_allsites_flow_water_mean <- readRDS(url("https://www.dropbox.com/scl/fi/n2zrmo8zmlo8rsvlv493e/fit_allsites_flow_water_mean.rds?rlkey=faftd5gfxrvgmmk88bg0doelv&dl=1"))
fit_allsites_flow_water_cv <- readRDS(url("https://www.dropbox.com/scl/fi/8c3xa0sovgrpolngdslqf/fit_allsites_flow_water_cv.rds?rlkey=jjx9v5izzhyzcuuvqtwgtauhw&dl=1"))

fit_allsites_spi_temp_mean <- readRDS(url("https://www.dropbox.com/scl/fi/frx00veu9q06wzmxzitbw/fit_allsites_spi_temp_mean.rds?rlkey=9b2efew2siyxwa6wm88c4uwxe&dl=1"))
fit_allsites_spi_temp_cv <- readRDS(url("https://www.dropbox.com/scl/fi/vc0saw1n88gyrkl2uwu4j/fit_allsites_spi_temp_cv.rds?rlkey=xuqq69qcgxwunha8yld2lbxl4&dl=1"))
fit_allsites_spi_water_mean <- readRDS(url("https://www.dropbox.com/scl/fi/yfblw00mhlpx2bn11l655/fit_allsites_spi_water_mean.rds?rlkey=a257vf3lv4uf1y42jlxqa86se&dl=1"))
fit_allsites_spi_water_cv <- readRDS(url("https://www.dropbox.com/scl/fi/9c1g64vktpipekcsyr8l0/fit_allsites_spi_water_cv.rds?rlkey=8t0zy6tjbldjmr2zwv8zmoazi&dl=1"))

fit_allsites_grow_temp_mean <- readRDS(url("https://www.dropbox.com/scl/fi/1i3e3yt01df9t7zwkzbew/fit_allsites_grow_temp_mean.rds?rlkey=7xy0kzfycilta76t7ixwjfune&dl=1"))
fit_allsites_grow_temp_cv <- readRDS(url("https://www.dropbox.com/scl/fi/vc0saw1n88gyrkl2uwu4j/fit_allsites_spi_temp_cv.rds?rlkey=xuqq69qcgxwunha8yld2lbxl4&dl=1"))
fit_allsites_grow_water_mean <- readRDS(url("https://www.dropbox.com/scl/fi/a2dprxvwrvzshv0xn45ww/fit_allsites_grow_water_mean.rds?rlkey=vt74d5ixqjzj4nbpllco1juxp&dl=1"))
fit_allsites_grow_water_cv <- readRDS(url("https://www.dropbox.com/scl/fi/jccvkcexr1e3mkjl6yrpx/fit_allsites_grow_water_cv.rds?rlkey=4cc0t8d5joi8aoqd3dbxk3rbh&dl=1"))

## Chains convergence
mcmc_trace(fit_allsites_surv_temp_mean, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                                         bendotemp_s,bherbtemp_s,bendoherb_s))+theme_bw()
mcmc_trace(fit_allsites_surv_temp_cv, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                                       bendotemp_s,bherbtemp_s,bendoherb_s))+theme_bw()
mcmc_trace(fit_allsites_surv_temp_mean, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                                         bendotemp_s,bherbtemp_s,bendoherb_s))+theme_bw()
mcmc_trace(fit_allsites_surv_water_mean,pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                                          bendotemp_s,bherbtemp_s,bendoherb_s))+theme_bw()
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
mcmc_trace(fit_allsites_spi_temp_cv, pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,btemp_spk,
                                                                      bendotemp_spk,bherbtemp_spk,bendoherb_spk))+theme_bw()
mcmc_trace(fit_allsites_spi_water_mean,pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,btemp_spk,
                                                                         bendotemp_spk,bherbtemp_spk,bendoherb_spk))+theme_bw()
mcmc_trace(fit_allsites_spi_water_cv, pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,btemp_spk,
                                                                       bendotemp_spk,bherbtemp_spk,bendoherb_spk))+theme_bw()
mcmc_trace(fit_allsites_grow_temp_mean, pars = quote_bare(b0_g,bendo_g,bherb_g,btemp_g,btemp_g,bendotemp_g,
                                                                         bherbtemp_g,bendoherb_g))+theme_bw()
mcmc_trace(fit_allsites_grow_temp_cv, pars = quote_bare(b0_g,bendo_g,bherb_g,btemp_g,btemp_g,bendotemp_g,
                                                                       bherbtemp_g,bendoherb_g))+theme_bw()
mcmc_trace(fit_allsites_grow_water_mean, pars = quote_bare(b0_g,bendo_g,bherb_g,btemp_g,btemp_g,bendotemp_g,
                                                                          bherbtemp_g,bendoherb_g))+theme_bw()
mcmc_trace(fit_allsites_grow_water_cv, pars = quote_bare(b0_g,bendo_g,bherb_g,btemp_g,btemp_g,bendotemp_g,
                                                                        bherbtemp_g,bendoherb_g))+theme_bw()


## Posterior predictive check for survival

## Model selection -----
## Model comparison based on epdl/looic 
### Survival
log_lik_surv_temp_mean <- loo::extract_log_lik(fit_allsites_surv_temp_mean, merge_chains = FALSE)
r_eff_surv_temp_mean <- loo::relative_eff(exp(log_lik_surv_temp_mean))
loo_surv_temp_mean <- loo(log_lik_surv_temp_mean, r_eff = r_eff_surv_temp_mean, cores = 1)
plot(loo_surv_temp_mean)

log_lik_surv_temp_cv <- loo::extract_log_lik(fit_allsites_surv_temp_cv, merge_chains = FALSE)
r_eff_surv_temp_cv <- loo::relative_eff(exp(log_lik_surv_temp_cv))
loo_surv_temp_cv <- loo(log_lik_surv_temp_cv, r_eff = r_eff_surv_temp_cv, cores = 3)
#plot(loo_surv_temp_cv)

log_lik_surv_water_mean <- loo::extract_log_lik(fit_allsites_surv_water_mean, merge_chains = FALSE)
r_eff_surv_water_mean <- loo::relative_eff(exp(log_lik_surv_water_mean))
loo_surv_water_mean <- loo(log_lik_surv_water_mean, r_eff = r_eff_surv_water_mean, cores = 3)
#plot(loo_surv_water_mean)

log_lik_surv_water_cv <- loo::extract_log_lik(fit_allsites_surv_water_cv, merge_chains = FALSE)
r_eff_surv_water_cv <- loo::relative_eff(exp(log_lik_surv_water_cv))
loo_surv_water_cv <- loo(log_lik_surv_water_cv, r_eff = r_eff_surv_water_cv, cores = 3)
#plot(loo_surv_water_cv)

(comp_surv <- loo::loo_compare(loo_surv_temp_mean,loo_surv_temp_cv, loo_surv_water_mean,loo_surv_water_cv))

### Growth 
log_lik_grow_temp_mean <- loo::extract_log_lik(fit_allsites_grow_temp_mean, merge_chains = FALSE)
r_eff_grow_temp_mean <- loo::relative_eff(exp(log_lik_grow_temp_mean))
loo_grow_temp_mean <- loo(log_lik_grow_temp_mean, r_eff = r_eff_grow_temp_mean, cores = 1)
# plot(loo_grow_temp_mean)

log_lik_grow_temp_cv <- loo::extract_log_lik(fit_allsites_grow_temp_cv, merge_chains = FALSE)
r_eff_grow_temp_cv <- loo::relative_eff(exp(log_lik_grow_temp_cv))
loo_grow_temp_cv <- loo(log_lik_grow_temp_cv, r_eff = r_eff_grow_temp_cv, cores = 3)
# plot(loo_grow_temp_cv)

log_lik_grow_water_mean <- loo::extract_log_lik(fit_allsites_grow_water_mean, merge_chains = FALSE)
r_eff_grow_water_mean <- loo::relative_eff(exp(log_lik_grow_water_mean))
loo_grow_water_mean <- loo(log_lik_grow_water_mean, r_eff = r_eff_grow_water_mean, cores = 3)
# plot(loo_grow_water_mean)

log_lik_grow_water_cv <- loo::extract_log_lik(fit_allsites_grow_water_cv, merge_chains = FALSE)
r_eff_grow_water_cv <- loo::relative_eff(exp(log_lik_grow_water_cv))
loo_grow_water_cv <- loo(log_lik_grow_water_cv, r_eff = r_eff_grow_water_cv, cores = 3)
# plot(loo_grow_water_cv)

(comp_grow <- loo::loo_compare(loo_grow_temp_mean,loo_grow_temp_cv, loo_grow_water_mean,loo_grow_water_cv))

## Flowering 
log_lik_flow_temp_mean <- loo::extract_log_lik(fit_allsites_flow_temp_mean, merge_chains = FALSE)
r_eff_flow_temp_mean <- loo::relative_eff(exp(log_lik_flow_temp_mean))
loo_flow_temp_mean <- loo(log_lik_flow_temp_mean, r_eff = r_eff_flow_temp_mean, cores = 1)
# plot(loo_flow_temp_mean)

log_lik_flow_temp_cv <- loo::extract_log_lik(fit_allsites_flow_temp_cv, merge_chains = FALSE)
r_eff_flow_temp_cv <- loo::relative_eff(exp(log_lik_flow_temp_cv))
loo_flow_temp_cv <- loo(log_lik_flow_temp_cv, r_eff = r_eff_flow_temp_cv, cores = 3)
# plot(loo_flow_temp_cv)

log_lik_flow_water_mean <- loo::extract_log_lik(fit_allsites_flow_water_mean, merge_chains = FALSE)
r_eff_flow_water_mean <- loo::relative_eff(exp(log_lik_flow_water_mean))
loo_flow_water_mean <- loo(log_lik_flow_water_mean, r_eff = r_eff_flow_water_mean, cores = 3)
# plot(loo_flow_water_mean)

log_lik_flow_water_cv <- loo::extract_log_lik(fit_allsites_flow_water_cv, merge_chains = FALSE)
r_eff_flow_water_cv <- loo::relative_eff(exp(log_lik_flow_water_cv))
loo_flow_water_cv <- loo(log_lik_flow_water_cv, r_eff = r_eff_flow_water_cv, cores = 3)
# plot(loo_flow_water_cv)

(comp_flow <- loo::loo_compare(loo_flow_temp_mean,loo_flow_temp_cv, loo_flow_water_mean,loo_flow_water_cv))

## Spiekelet
log_lik_spi_temp_mean <- loo::extract_log_lik(fit_allsites_spi_temp_mean, merge_chains = FALSE)
r_eff_spi_temp_mean <- loo::relative_eff(exp(log_lik_spi_temp_mean))
loo_spi_temp_mean <- loo(log_lik_spi_temp_mean, r_eff = r_eff_spi_temp_mean, cores = 1)
#plot(loo_spi_temp_mean)

log_lik_spi_temp_cv <- loo::extract_log_lik(fit_allsites_spi_temp_cv, merge_chains = FALSE)
r_eff_spi_temp_cv <- loo::relative_eff(exp(log_lik_spi_temp_cv))
loo_spi_temp_cv <- loo(log_lik_spi_temp_cv, r_eff = r_eff_spi_temp_cv, cores = 3)
#plot(loo_spi_temp_cv)

log_lik_spi_water_mean <- loo::extract_log_lik(fit_allsites_spi_water_mean, merge_chains = FALSE)
r_eff_spi_water_mean <- loo::relative_eff(exp(log_lik_spi_water_mean))
loo_spi_water_mean <- loo(log_lik_spi_water_mean, r_eff = r_eff_spi_water_mean, cores = 3)
#plot(loo_spi_water_mean)

log_lik_spi_water_cv <- loo::extract_log_lik(fit_allsites_spi_water_cv, merge_chains = FALSE)
r_eff_spi_water_cv <- loo::relative_eff(exp(log_lik_spi_water_cv))
loo_spi_water_cv <- loo(log_lik_spi_water_cv, r_eff = r_eff_spi_water_cv, cores = 3)
#plot(loo_spi_water_cv)

(comp_spi <- loo::loo_compare(loo_spi_temp_mean,loo_spi_temp_cv, loo_spi_water_mean,loo_spi_water_cv))

## # Posterior mean values for each vital rate----
posterior_surv <- as.array(fit_allsites_surv_temp_cv)
color_scheme_set("red")
surv<-mcmc_intervals(posterior_surv, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,
                                                       bendotemp_s,bherbtemp_s,bendoherb_s)) + 
  ggplot2::scale_y_discrete(limits = c("b0_s","bendo_s","bherb_s","btemp_s",
                                       "bendotemp_s","bherbtemp_s","bendoherb_s"),
                            labels=c("b0_s"="Intercept",
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
                            labels=c("b0_g"="Intercept",
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
                            labels=c("b0_f"="Intercept",
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

posterior_spi <- as.array(fit_allsites_spi_temp_cv)
flow<-mcmc_intervals(posterior_spi, pars = quote_bare(b0_sp,bendo_f,bherb_f,btemp_f,
                                                       bendotemp_f,bherbtemp_f,bendoherb_f)) + 
  ggplot2::scale_y_discrete(limits = c("b0_f","bendo_f","bherb_f","btemp_f",
                                       "bendotemp_f","bherbtemp_f","bendoherb_f"),
                            labels=c("b0_f"="Intercept",
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
  xlim(-5,5)+
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

posterior_spi <- as.array(fit_allsites_spi_temp_cv)
spi<-mcmc_intervals(posterior_spi, pars = quote_bare(b0_spk,bendo_spk,bherb_spk,btemp_spk,
                                                       bendotemp_spk,bherbtemp_spk,bendoherb_spk)) + 
  ggplot2::scale_y_discrete(limits = c("b0_spk","bendo_spk","bherb_spk","btemp_spk",
                                       "bendotemp_spk","bherbtemp_spk","bendoherb_spk"),
                            labels=c("b0_spk"="Intercept",
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
