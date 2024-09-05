## Author: Jacob Moutouama and Tom Miller 
## Purpose: Fit vital rate models to test the effect of grass-endophyte symbiosis and endophyte hyphal density on population demography 

rm(list = ls())

# Load required package
library(tidyverse)
library(readxl) 
library(bayesplot)
library(ggsci)
library(rstan)
# set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(13)

# quote a series of bare names
quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

invlogit<-function(x){exp(x)/(1+exp(x))}

# tom_path<-"C:/Users/tm9/Dropbox/github/ELVI-endophyte-density" 
jacob_path<-"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density"

# HOBO data ----
## format date and separate year-month-day
choose_path<-jacob_path
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

# pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/climatesite.pdf",height =5,width =12,useDingbats = F)
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

## how strongly are these summary stats related to duration of the experimental period?
par(mfrow=c(2,2))
plot(HOBO_summary$duration,HOBO_summary$temp_mean,pch=16)
plot(HOBO_summary$duration,HOBO_summary$temp_cv,pch=16)
plot(HOBO_summary$duration,HOBO_summary$water_mean,pch=16)
plot(HOBO_summary$duration,HOBO_summary$water_cv,pch=16)


# summary(HOBO_summary)
#unique(HOBO_summary$Site)
# Correlation between temperature and soil moisture

(figcorrelation_mean<-ggplot(data=HOBO_summary, aes(x=scale(temp_mean), y=scale(water_mean))) +
  geom_smooth(method="lm",color = "blue") +
  geom_point(color = 'blue',shape = 19, size=4) +
  labs (x="Water (mean)",y="Temperature (mean)") +
  ggpubr::stat_cor(aes(label=..rr.label..)) +
  theme_bw())

(figcorrelation_cv<-ggplot(data=HOBO_summary, aes(x=scale(temp_cv), y=scale(water_cv))) +
  geom_smooth(method="lm",color = "blue") +
  geom_point(color = 'blue',shape = 19, size=4) +
  labs (x="Water (CV)",y="Temperature (CV)") +
  ggpubr::stat_cor(aes(label=..rr.label..)) +
  theme_bw())

# pdf("/Users/jmoutouama/Dropbox/Miller Lab/github/ELVI-endo-density/Figure/correlation.pdf",height=3,width =6,useDingbats = F)
# (figcorrelation<-ggpubr::ggarrange(figcorrelation_mean,figcorrelation_cv,labels = c("A","B")))
# dev.off()

# barplot(HOBO_spring_summary$water_mean,
#         names.arg = HOBO_spring_summary$Site)
# barplot(HOBO_spring_summary$temp_mean,
#         names.arg = HOBO_spring_summary$Site)

# Demographic data -----
# Merge the demographic with the climatic data
datini<-read_csv(paste0(choose_path,"/Data/Initialdata.csv"))
dat23<-read_csv(paste0(choose_path,"/Data/census2023.csv"))
dat24<-read_csv(paste0(choose_path,"/Data/census2024.csv"))
datherbivory<-read_csv(paste0(choose_path,"/Data/herbivory.csv"))

# calculate the total spikelet for each census
dat23 %>% 
  mutate(Spikelet_23=rowSums(across(Spikelet_A:Spikelet_C)))->dat23_spike
dat24 %>% 
  mutate(Spikelet_24=rowSums(across(Spikelet_A:Spikelet_C)),Inf_24=rowSums(across(attachedInf_24:brokenInf_24)))->dat24_spike_inf

## Merge the initial data with the 23 data and the 23 data with the 24 -----
datini23 <- left_join(x = datini,y =dat23_spike,by=c("Tag_ID"))
## warning says this id in x occurs multiple times in y
datini$Tag_ID[1455]
datini$Tag_ID[628]
dat23[dat23$Tag_ID==datini$Tag_ID[1455],] ##this needs to be fixed
dat23[dat23$Tag_ID==datini$Tag_ID[628],] ##this needs to be fixed

datini23 %>% 
  mutate(tiller_t=ini_Tiller,
         tiller_t1=Tiller_23,
         inf_t1=Inf_23,
         spikelet_t1=Spikelet_23,
         year=c(rep(2023,nrow(datini23)))) %>% 
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
                tiller_Herb,
                year)->datini23_t_t1

dat2324 <- left_join(x = datini23 ,y =dat24_spike_inf,by=c("Tag_ID")) 

dat2324 %>% 
  mutate(tiller_t=Tiller_23,
         tiller_t1=Tiller_24,
         inf_t1=Inf_24,
         spikelet_t1=Spikelet_24,
         tiller_Herb=tiller_Herb.y,
         year=c(rep(2024,nrow(dat2324)))) %>% 
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
                tiller_Herb,
                year)->dat2324_t_t1

datdemo<-rbind(datini23_t_t1,dat2324_t_t1)
# names(datdemo)
# unique(datdemo$Species)
## Merge the demographic data with the herbivory data -----
demography<-left_join(x=datdemo,y=datherbivory,by=c("Site","Plot","Species"))# Merge the demographic data with the herbivory data
# view(demography)
# head(demography)
HOBO_summary %>% 
  rename(Site=site)->HOBO_summary_clean
## Merge the demographic data with the climatic data -----
demography_climate<-left_join(x=demography,y=HOBO_summary_clean,by=c("Site"))# Merge the demographic data with the temperature data
# Subset only ELVI data -----
demography_climate %>% 
  filter(Species=="ELVI")->demography_climate_elvi
demography_climate_elvi$surv1<-1*(!is.na(demography_climate_elvi$tiller_t) & !is.na(demography_climate_elvi$tiller_t1))
demography_climate_elvi$site_plot<-interaction(demography_climate_elvi$Site,demography_climate_elvi$Plot)
# names(demography_climate_elvi)
# view(demography_climate_elvi)
# summary(demography_climate_elvi)
# H1: We hypothesized that stress associated with aridity and low precipitation would strengthen the plant-fungal mutualism, such that the fitness benefits of endophyte symbiosis are maximized at the range edge. 

## Survival 
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

## Chains convergence
traceplot(fit_allsites_surv_temp_mean,inc_warmup=TRUE, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                 bendotemp_s,bherbtemp_s,bendoherb_s))+theme_bw()
pairs(fit_allsites_surv_temp_mean, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                           bendotemp_s,bherbtemp_s,bendoherb_s),las=1)
traceplot(fit_allsites_surv_temp_cv,inc_warmup=TRUE, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                         bendotemp_s,bherbtemp_s,bendoherb_s))+theme_bw()
pairs(fit_allsites_surv_temp_cv, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                     bendotemp_s,bherbtemp_s,bendoherb_s),las=1)
traceplot(fit_allsites_surv_temp_mean,inc_warmup=TRUE, pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                         bendotemp_s,bherbtemp_s,bendoherb_s))+theme_bw()
pairs(fit_allsites_surv_temp_mean,pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                     bendotemp_s,bherbtemp_s,bendoherb_s),las=1)
traceplot(fit_allsites_surv_water_mean, inc_warmup=TRUE,pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                         bendotemp_s,bherbtemp_s,bendoherb_s))+theme_bw()
pairs(fit_allsites_surv_water_mean,pars = quote_bare(b0_s,bendo_s,bherb_s,btemp_s,btemp_s,
                                                     bendotemp_s,bherbtemp_s,bendoherb_s),las=1)
## Save RDS file for further use
# saveRDS(fit_allsites_surv_temp_mean, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/output/Survival/fit_allsites_surv_temp_mean.rds')
# saveRDS(fit_allsites_surv_temp_cv, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/output/Survival/fit_allsites_surv_temp_cv.rds')
# saveRDS(fit_allsites_surv_water_mean, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/output/Survival/fit_allsites_surv_water_mean.rds')
# saveRDS(fit_allsites_surv_water_cv, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/output/Survival/fit_allsites_surv_water_cv.rds')

#load stan output 
fit_allsites_surv_temp_mean <- readRDS(url("https://www.dropbox.com/scl/fi/g3spz9lff9xy3153o96vo/fit_allsites_surv_temp_mean.rds?rlkey=jnv64cw87w1wi9te8brvfh7df&dl=1"))
fit_allsites_surv_temp_cv <- readRDS(url("https://www.dropbox.com/scl/fi/vq8ryu7vgtp3yqs68ewff/fit_allsites_surv_temp_cv.rds?rlkey=i227tu61cp224hkjpgekox11r&dl=1"))
fit_allsites_surv_water_mean <- readRDS(url("https://www.dropbox.com/scl/fi/at18rlo9e53pqc1tmmund/fit_allsites_surv_water_mean.rds?rlkey=hvvbxzn9g1x4brlidxcm147n3&dl=1"))
fit_allsites_surv_water_cv <- readRDS(url("https://www.dropbox.com/scl/fi/wdnd6e2oc4zunmmhsvayn/fit_allsites_surv_water_cv.rds?rlkey=sbim89fp48gnsoagpk1r57t4k&dl=1"))

## Model comparison based on epdl/looic 
log_lik_surv_temp_mean <- loo::extract_log_lik(fit_allsites_surv_temp_mean, merge_chains = FALSE)
r_eff_surv_temp_mean <- loo::relative_eff(exp(log_lik_surv_temp_mean))
loo_surv_temp_mean <- loo(log_lik_surv_temp_mean, r_eff = r_eff_surv_temp_mean, cores = 1)
plot(loo_surv_temp_mean)

log_lik_surv_temp_cv <- loo::extract_log_lik(fit_allsites_surv_temp_cv, merge_chains = FALSE)
r_eff_surv_temp_cv <- loo::relative_eff(exp(log_lik_surv_temp_cv))
loo_surv_temp_cv <- loo(log_lik_surv_temp_cv, r_eff = r_eff_surv_temp_cv, cores = 3)
plot(loo_surv_temp_cv)

log_lik_surv_water_mean <- loo::extract_log_lik(fit_allsites_surv_water_mean, merge_chains = FALSE)
r_eff_surv_water_mean <- loo::relative_eff(exp(log_lik_surv_water_mean))
loo_surv_water_mean <- loo(log_lik_surv_water_mean, r_eff = r_eff_surv_water_mean, cores = 3)
plot(loo_surv_water_mean)

log_lik_surv_water_cv <- loo::extract_log_lik(fit_allsites_surv_water_cv, merge_chains = FALSE)
r_eff_surv_water_cv <- loo::relative_eff(exp(log_lik_surv_water_cv))
loo_surv_water_cv <- loo(log_lik_surv_water_cv, r_eff = r_eff_surv_water_cv, cores = 3)
plot(loo_surv_water_cv)

(comp_surv <- loo::loo_compare(loo_surv_temp_mean,loo_surv_temp_cv, loo_surv_water_mean,loo_surv_water_cv))

## Flowering 
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




traceplot(fit_allsites_flow_temp_mean,inc_warmup=TRUE, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,btemp_f,
                                                                         bendotemp_f,bherbtemp_f,bendoherb_f))+theme_bw()
pairs(fit_allsites_flow_temp_mean, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,btemp_f,
                                                     bendotemp_f,bherbtemp_f,bendoherb_f),las=1)
traceplot(fit_allsites_flow_temp_cv,inc_warmup=TRUE, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,btemp_f,
                                                                         bendotemp_f,bherbtemp_f,bendoherb_f))+theme_bw()
pairs(fit_allsites_flow_temp_cv, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,btemp_f,
                                                     bendotemp_f,bherbtemp_f,bendoherb_f),las=1)
traceplot(fit_allsites_flow_water_mean,inc_warmup=TRUE, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,btemp_f,
                                                                       bendotemp_f,bherbtemp_f,bendoherb_f))+theme_bw()
pairs(fit_allsites_flow_water_mean, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,btemp_f,
                                                   bendotemp_f,bherbtemp_f,bendoherb_f),las=1)
traceplot(fit_allsites_flow_water_cv,inc_warmup=TRUE, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,btemp_f,
                                                                          bendotemp_f,bherbtemp_f,bendoherb_f))+theme_bw()
pairs(fit_allsites_flow_water_cv, pars = quote_bare(b0_f,bendo_f,bherb_f,btemp_f,btemp_f,
                                                      bendotemp_f,bherbtemp_f,bendoherb_f),las=1)

# Save RDS file for further use
# saveRDS(fit_allsites_flow_temp_mean, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/output/Flowering/fit_allsites_flow_temp_mean.rds')
# saveRDS(fit_allsites_flow_temp_cv, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/output/Flowering/fit_allsites_flow_temp_cv.rds')
# saveRDS(fit_allsites_flow_water_mean, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/output/Flowering/fit_allsites_flow_water_mean.rds')
# saveRDS(fit_allsites_flow_water_cv, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/output/Flowering/fit_allsites_flow_water_cv.rds')


## Read and format the growth data to build the model
demography_climate_elvi$grow<-(log(demography_climate_elvi$tiller_t1+1) - log(demography_climate_elvi$tiller_t+1))# Relative growth rate
hist(demography_climate_elvi$grow)

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
summary(data_sites_grow_temp_mean$y_g)
data_sites_grow_temp_cv <- list( n_sites    = demography_climate_elvi_grow$Site %>% n_distinct,
                                 n_pops  = demography_climate_elvi_grow$Population %>% n_distinct(),
                                 # survival data
                                 n_plot_g = demography_climate_elvi_grow$Plot %>% n_distinct,
                                 site_s   = demography_climate_elvi_grow$Site,
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

fit_allsites_grow_temp_mean <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Analysis/stan/elvi_growth.stan",
  data = data_sites_grow_temp_mean,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)

traceplot(fit_allsites_grow_temp_mean,inc_warmup=TRUE, pars = quote_bare(b0_g,bendo_g,bherb_g,btemp_g,btemp_g,
                                                                         bendotemp_g,bherbtemp_g,bendoherb_g))+theme_bw()
pairs(fit_allsites_grow_temp_mean, pars = quote_bare(b0_g,bendo_g,bherb_g,btemp_g,btemp_g,
                                                     bendotemp_g,bherbtemp_g,bendoherb_g),las=1)


demography_climate_elvi %>% 
  subset( tiller_t > 0 & tiller_t1 > 0)%>%
  dplyr::select(year, Population, Site, Plot,site_plot, Endo, Herbivory,
                tiller_t, grow,temp_mean,temp_cv,water_mean,water_cv)%>% 
  na.omit %>% 
  mutate( Site= Site %>% as.factor,
          Plot = Plot %>% as.factor,
          site_plot=site_plot %>% as.factor ,
          Endo = Endo %>% as.factor,
          Herbivory=Herbivory %>% as.factor,
          Population = Population %>% as.factor) %>%
  mutate( log_size_t0 = log(tiller_t),
          grow_t1=grow,
          log_temp_mean = log(temp_mean),
          log_temp_cv = log(temp_cv),
          log_water_mean = log(water_mean),
          log_water_cv = log(water_cv))->demography_climate_elvi_grow_brms


modg = brms::brm(
  grow_t1 ~ log_temp_mean + Endo + Herbivory + log_temp_mean*Endo + log_temp_mean*Herbivory + (1|Population) + (1|Site/Plot), data = demography_climate_elvi_grow_brms,
  family = 'gaussian',
  prior = set_prior('normal(0, 3)'), iter = 1000,
  chains = 4,
  cores = 4
)

plot(modg)
