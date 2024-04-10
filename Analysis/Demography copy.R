## Author: Jacob Moutouama and Tom Miller 
## Purpose: Fit statistical models to vital rate to test the effect of demographic effects of grass-endophyte symbiosis and endophyte hyphal density 

rm(list = ls())

# Load required package

library(tidyverse)
library(readxl) 
library(lme4)
library(lmerTest)
library(bbmle)
library(glmmTMB)
library(ggsci)

# tools::package_dependencies("Matrix", which = "LinkingTo", reverse = TRUE)[[1L]]
# install.packages("lme4", type = "source")
# oo <- options(repos = "https://cran.r-project.org/")
# install.packages("Matrix")
# install.packages("lme4")
# options(oo)
# choose path for data import

# tom_path<-"G:/Shared drives/Miller Lab/Endophytes - Range Limits/Elymus Data Analysis" # You need to change this because I have everything on Github
jacob_path<-"/Users/jmoutouama/Dropbox/Miller Lab/github/ELVI-endo-density"

# HOBO data ----
## format date and separate year-month-day
choose_path<-jacob_path
list.files(path = paste0(choose_path,"/Data/HOBO data/"),  
           pattern = "*.xlsx", full.names = TRUE) %>% # Identify all excel files
  lapply(read_excel) %>%                              # Store all files in list
  bind_rows ->hobo_data_raw # get HOBO data

tidyr::separate(hobo_data_raw, 'date',
                into = c('longdate', 'time'),
                sep= ' ') %>%
  tidyr::separate('longdate', # Separate the ‘longdate’ column into separate columns for month, day and year using the separate() function.
                  into = c('year','month', 'day'),
                  sep= '-',
                  remove = FALSE)->hobo_data_full 

## double check if the starting and ending dates are correct
hobo_data_full %>% 
  group_by(Site) %>% 
  summarise(start=range(longdate)[1],
            end=range(longdate)[2]) 

## average over days to look at overall trend across sites
hobo_data_full %>% 
  group_by(longdate,Site,day) %>% 
  summarise(daily_mean_moist=mean(water),daily_mean_temp=mean(temp))->HOBO_daily


## Plot the daily trend for temperature and soil moisture from start to end

figtempsite<-ggplot(HOBO_daily, aes(x=as.Date(longdate, format= "%Y - %m - %d"), y=daily_mean_temp))+
  geom_line(aes(colour=Site))+
  ggtitle("a")+
  scale_fill_jco()+
  theme_bw()+ 
  labs( y="Daily temperature  (°C)", x="Month") 

figmoistsite<-ggplot(HOBO_daily, aes(x=as.Date(longdate, format= "%Y - %m - %d"), y=daily_mean_moist))+
  geom_line(aes(colour=Site))+
  ggtitle("b")+
  scale_fill_jco()+
  theme_bw()+ 
  labs( y="Daily soil moisture (wfv)", x="Month") # Plot the raw trend from start to end

# pdf("/Users/jmoutouama/Dropbox/Miller Lab/github/ELVI-endo-density/Figure/climatesite.pdf",height =3.5,width =6,useDingbats = F)
# (Figclimatesite<-ggpubr::ggarrange(figtempsite,figmoistsite,common.legend = TRUE,legend.position = "bottom",legend.box = "horizontal"))
# dev.off()

## Average over days

HOBO_daily %>%
  group_by(Site)%>%
  summarise(water_mean = mean(daily_mean_moist),
            water_cv=sd(daily_mean_moist)/water_mean,
            temp_mean = mean(daily_mean_temp),
            temp_cv=sd(daily_mean_temp)/temp_mean)->HOBO_summary
summary(HOBO_summary)
#unique(HOBO_summary$Site)
# Correlation between temperature and soil moisture

(figcorrelation_mean<-ggplot(data=HOBO_summary, aes(x=scale(temp_mean), y=scale(water_mean))) +
  geom_smooth(method="lm",color = "blue") +
  geom_point(color = 'blue',shape = 19, size=4) +
  labs (x="Water (mean)",y="Temperature (mean)") +
  ggpubr::stat_cor(aes(label=..rr.label..)) +
  theme_bw())
dev.off()

figcorrelation_cv<-ggplot(data=HOBO_summary, aes(x=scale(temp_cv), y=scale(water_cv))) +
  geom_smooth(method="lm",color = "blue") +
  geom_point(color = 'blue',shape = 19, size=4) +
  labs (x="Water (CV)",y="Temperature (CV)") +
  ggpubr::stat_cor(aes(label=..rr.label..)) +
  theme_bw()

pdf("/Users/jmoutouama/Dropbox/Miller Lab/github/ELVI-endo-density/Figure/correlation.pdf",height=3,width =6,useDingbats = F)
(figcorrelation<-ggpubr::ggarrange(figcorrelation_mean,figcorrelation_cv,labels = c("A","B")))
dev.off()

# barplot(HOBO_spring_summary$water_mean,
#         names.arg = HOBO_spring_summary$Site)
# barplot(HOBO_spring_summary$temp_mean,
#         names.arg = HOBO_spring_summary$Site)

# Demographic data -----
# Merge the demographic with the climatic data
datini<-read_csv(paste0(choose_path,"/Data/Initialdata.csv"))
dat23<-read_csv(paste0(choose_path,"/Data/census2023.csv"))
# datherbivory<-read_csv(paste0(choose_path,"/Data/herbivory.csv"))
# head(datherbivory)
## Merge the initial data with the 23 data -----
datdemo <- left_join(x = datini ,y =dat23,by=c("Tag_ID")) 
# names(dat)
# unique(dat$Species)
#demography<-left_join(x=datdemo,y=datherbivory,by=c("Site","Species","Plot"))# Merge the demographic data with the herbivory data
# head(demography)
demography_climate<-left_join(x=datdemo,y=HOBO_summary,by=c("Site"))# Merge the demographic data with the temperature data
demography_climate %>% 
  filter(Species=="ELVI")->demography_climate_elvi
#view(demography_climate_elvi)
# H1: We hypothesized that stress associated with aridity and low precipitation would strengthen the plant-fungal mutualism, such that the fitness benefits of endophyte symbiosis are maximized at the range edge. 

demography_climate_elvi %>% 
  filter(Tiller_23 >0)->demography_climate_elvi
demography_climate_elvi$grow<-(log(demography_climate_elvi$Tiller_23) - log(demography_climate_elvi$ini_Tiller))# Relative growth rate
#summary(demography_tempelvi)
demography_climate_elvi$Endo<-as.factor(demography_climate_elvi$Endo)
demography_climate_elvi %>% 
  mutate(genotype = interaction(Population,GreenhouseID))->demography_climate_elvi
# table(demography_tempelvi$genotype,demography_tempelvi$Site)
view(demography_climate_elvi)
# Soil moisture
## Growth----

endostatus_growth_soil_mods <- list()
endostatus_growth_soil_mods[[1]] <-lmer(grow ~ 1 + (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_soil_mods[[2]] <- lmer(grow ~  Endo + (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_soil_mods[[3]] <- lmer(grow ~ Endo + water_mean + (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_soil_mods[[4]] <- lmer(grow ~ Endo*water_mean + (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_soil_mods[[5]] <- lmer(grow ~ Endo + water_cv+ (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_soil_mods[[6]] <- lmer(grow ~ Endo*water_cv+ (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)

AICtab(endostatus_growth_soil_mods)
summary(endostatus_growth_soil_mods[[5]])

# ranef(endostatus_growth_mods[[5]])$Site

coefsoilgrowth<-fixef(endostatus_growth_soil_mods[[5]])
s.igEndo0<-coefsoilgrowth[1]			## intercept for Endo 0
s.sgEndo0<-coefsoilgrowth[3]			## slope for for Endo 0

s.igEndo1 <-coefsoilgrowth[2]+s.igEndo0	## intercept for Endo 1
s.sgEndo1 <-coefsoilgrowth[4]+s.sgEndo0	## slope for Endo 1


# Number inflorescence ---- 
endostatus_inf_soil_mods<-list()
endostatus_inf_soil_mods[[1]] <- glmmTMB(Inf_23 ~ 1+ (1|Site/Plot) + (1|Population) ,ziformula = ~1,family = nbinom2, data=demography_tempelvi,na.action = na.omit)
endostatus_inf_soil_mods[[2]] <- glmmTMB(Inf_23 ~ Endo + (1|Site/Plot) + (1|Population),ziformula = ~1,family = nbinom2,  data=demography_tempelvi)
endostatus_inf_soil_mods[[3]] <- glmmTMB(Inf_23 ~ Endo + water_mean + (1|Site/Plot) + (1|Population) ,ziformula = ~1,family = nbinom2,  data=demography_tempelvi)
endostatus_inf_soil_mods[[4]] <- glmmTMB(Inf_23 ~ Endo*water_mean + (1|Site/Plot) + (1|Population) ,ziformula = ~1,family = nbinom2,  data=demography_tempelvi)
endostatus_inf_soil_mods[[5]] <- glmmTMB(Inf_23 ~ Endo + water_cv + (1|Site/Plot) + (1|Population), ziformula = ~1,family = nbinom2, data=demography_tempelvi)
endostatus_inf_soil_mods[[6]] <- glmmTMB(Inf_23 ~ Endo*water_cv + (1|Site/Plot) + (1|Population) ,ziformula = ~1,family = nbinom2,  data=demography_tempelvi)

AICtab(endostatus_inf_soil_mods)
summary(endostatus_inf_soil_mods[[4]])

# coefinf<-fixef(endostatus_inf_mods[[4]])
# s.iiEndo0<-coefinf[1]$cond[1]			## intercept for Endo 0
# s.siEndo0<-coefinf[1]$cond[3]			## slope for for Endo 0
# 
# s.iiEndo1 <-coefinf[1]$cond[2]+ s.iiEndo0	## intercept for Endo 1
# s.siEndo1 <- s.siEndo0	## slope for Endo 1

# Number of spikelet 
endostatus_spike_soil_mods<-list()
endostatus_spike_soil_mods[[1]] <- glmmTMB(Spikelet ~ 1 + (1|Site/Plot) + (1|Population)  , data=demography_tempelvi,ziformula = ~1,family = nbinom2)
endostatus_spike_soil_mods[[2]] <- glmmTMB(Spikelet ~ Endo + (1|Site/Plot) + (1|Plot/Population)  ,ziformula = ~1,family = nbinom2, data=demography_tempelvi)
endostatus_spike_soil_mods[[3]] <- glmmTMB(Spikelet ~ Endo + water_mean + (1|Site/Plot) + (1|Plot/Population) ,ziformula = ~1,family = nbinom2, data=demography_tempelvi,control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
endostatus_spike_soil_mods[[4]] <- glmmTMB(Spikelet ~ Endo*water_mean + (1|Site/Plot) + (1|Plot/Population) ,ziformula = ~1,family = nbinom2, data=demography_tempelvi,control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
endostatus_spike_soil_mods[[5]] <- glmmTMB(Spikelet ~ Endo + water_cv + (1|Site/Plot) + (1|Plot/Population)  ,ziformula = ~1,family = nbinom2, data=demography_tempelvi)
endostatus_spike_soil_mods[[6]] <- glmmTMB(Spikelet ~ Endo*water_cv + (1|Site/Plot) + (1|Plot/Population) , ziformula = ~1,family = nbinom2,data=demography_tempelvi)
AICtab(endostatus_spike_soil_mods)
summary(endostatus_spike_soil_mods[[3]])



# coefspk<-fixef(endostatus_spike_mods[[3]])
# s.isEndo0<-coefspk[1]$cond[1]			## intercept for Endo 0
# s.ssEndo0<-coefspk[1]$cond[3]			## slope for for Endo 0
# 
# s.isEndo1 <-coefspk[1]$cond[2]+ s.isEndo0	## intercept for Endo 1
# s.ssEndo1 <- s.ssEndo0	## slope for Endo 1
# 


# Plot the relationships between each vital rate and the water CV/ Water mean 
cbPalette <- c("#999999", "#E69F00")
water_cv_seq<-seq(min(na.omit(demography_tempelvi$water_cv)),max(na.omit(demography_tempelvi$water_cv)),length.out=20)
water_mean_seq<-seq(min(na.omit(demography_tempelvi$water_mean)),max(na.omit(demography_tempelvi$water_mean)),length.out=20)

# layout.matrix <- rbind(matrix(1:3, nrow = 1, ncol = 3, byrow = F))
#print figure
pdf("/Users/jmoutouama/Dropbox/Miller Lab/github/ELVI-endo-density/Figure/Soilmoisture.pdf",height=3,width =9,useDingbats = F)
# layout(mat = layout.matrix,
#        heights = rep(c(1, 1,1),2), # Heights of the two rows
#        widths = c(3, 3,3))
# layout.show(3)
# par(oma=c(5,1,0.5,0.5))

par(mfrow=c(1,3),mar=c(5,5,2,1))
with(demography_tempelvi,{
  # par(mar=c(4,4,2,4))
  plot(water_cv,grow,cex.lab=1.5,xlab="Soil water content (CV)",ylab="Growth rate");box()
  points(jitter(water_cv),jitter(grow),bg=cbPalette,pch=21,col=NA,cex=2,lwd=2)
  title("A",adj=0)
  lines(water_cv_seq,(s.igEndo0 + s.sgEndo0*water_cv_seq),col=cbPalette[1],lwd=3,)
  lines(water_cv_seq,(s.igEndo1 + s.sgEndo1*water_cv_seq),col=cbPalette[2],lwd=3,)
})

with(demography_tempelvi,{
  # par(mar=c(4,4,2,4))
  plot(water_cv,Inf_23,cex.lab=1.5,xlab="Soil water content (CV)",ylab="# Inflorescences");box()
  points(jitter(water_cv),jitter(Inf_23),bg=cbPalette,pch=21,col=NA,cex=2,lwd=2)
  title("B",adj=0)
  # lines(water_cv_seq,exp(s.iiEndo0 + s.siEndo0*water_cv_seq),col=cbPalette[1],lwd=3,)
  # lines(water_cv_seq,exp(s.iiEndo1 + s.siEndo1*water_cv_seq),col=cbPalette[2],lwd=3,)
})

with(demography_tempelvi,{
  # par(mar=c(4,4,2,4))
  plot(water_mean,Spikelet,cex.lab=1.5,xlab="Soil water content (mean)",ylab="# Spikelets");box()
  points(jitter(water_mean),jitter(Spikelet),bg=cbPalette,pch=21,col=NA,cex=2,lwd=2)
  title("C",adj=0)
  # lines(water_mean_seq,exp(s.isEndo0 + s.ssEndo0*water_mean_seq),col=cbPalette[1],lwd=3,)
  # lines(water_mean_seq,exp(s.isEndo1 + s.ssEndo1*water_mean_seq),col=cbPalette[2],lwd=3,)
  legend("topright",legend=c("E-","E+"),lwd=2,col =cbPalette ,bty="n",cex=1)
})

dev.off()

# Temperature
## Growth----
# intercept and water_cv sig
endostatus_growth_temp_mods <- list()
endostatus_growth_temp_mods[[1]] <-lmer(grow ~ 1 + (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_temp_mods[[2]] <- lmer(grow ~  Endo + (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_temp_mods[[3]] <- lmer(grow ~ Endo + temp_mean + (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_temp_mods[[4]] <- lmer(grow ~ Endo*temp_mean + (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_temp_mods[[5]] <- lmer(grow ~ Endo + temp_cv+ (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_temp_mods[[6]] <- lmer(grow ~ Endo*temp_cv+ (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)

AICtab(endostatus_growth_temp_mods)
summary(endostatus_growth_temp_mods[[6]])

# ranef(endostatus_growth_mods[[5]])$Site

coefgrowth<-fixef(endostatus_growth_temp_mods[[5]])
s.igEndo0<-coefgrowth[1]			## intercept for Endo 0
s.sgEndo0<-coefgrowth[3]			## slope for for Endo 0

s.igEndo1 <-coefgrowth[2]+s.igEndo0	## intercept for Endo 1
s.sgEndo1 <-coefgrowth[4]+s.sgEndo0	## slope for Endo 1


# Number inflorescence ---- 
endostatus_inf_temp_mods<-list()
endostatus_inf_temp_mods[[1]] <- glmmTMB(Inf_23 ~ 1+ (1|Site/Plot) + (1|Population) ,ziformula = ~1,family = nbinom2, data=demography_tempelvi,na.action = na.omit)
endostatus_inf_temp_mods[[2]] <- glmmTMB(Inf_23 ~ Endo + (1|Site/Plot) + (1|Population),ziformula = ~1,family = nbinom2,  data=demography_tempelvi)
endostatus_inf_temp_mods[[3]] <- glmmTMB(Inf_23 ~ Endo + temp_mean + (1|Site/Plot) + (1|Population) ,ziformula = ~1,family = nbinom2,  data=demography_tempelvi)
endostatus_inf_temp_mods[[4]] <- glmmTMB(Inf_23 ~ Endo + temp_cv + (1|Site/Plot) + (1|Population), ziformula = ~1,family = nbinom2, data=demography_tempelvi)
endostatus_inf_temp_mods[[5]] <- glmmTMB(Inf_23 ~ Endo*temp_cv + (1|Site/Plot) + (1|Population) ,ziformula = ~1,family = nbinom2,  data=demography_tempelvi)

AICtab(endostatus_inf_temp_mods)
summary(endostatus_inf_temp_mods[[4]])

coefinf<-fixef(endostatus_inf_temp_mods[[4]])
s.iiEndo0<-coefinf[1]$cond[1]			## intercept for Endo 0
s.siEndo0<-coefinf[1]$cond[3]			## slope for for Endo 0

s.iiEndo1 <-coefinf[1]$cond[2]+ s.iiEndo0	## intercept for Endo 1
s.siEndo1 <- s.siEndo0	## slope for Endo 1

# Number of spikelet 
endostatus_spike_temp_mods<-list()
endostatus_spike_temp_mods[[1]] <- glmmTMB(Spikelet ~ 1 + (1|Site/Plot) + (1|Population)  , data=demography_tempelvi,ziformula = ~1,family = nbinom2)
endostatus_spike_temp_mods[[2]] <- glmmTMB(Spikelet ~ Endo + (1|Site/Plot) + (1|Plot/Population)  ,ziformula = ~1,family = nbinom2, data=demography_tempelvi)
endostatus_spike_temp_mods[[3]] <- glmmTMB(Spikelet ~ Endo + temp_mean + (1|Site/Plot) + (1|Plot/Population) ,ziformula = ~1,family = nbinom2, data=demography_tempelvi,control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
endostatus_spike_temp_mods[[4]] <- glmmTMB(Spikelet ~ Endo + temp_cv + (1|Site/Plot) + (1|Plot/Population)  ,ziformula = ~1,family = nbinom2, data=demography_tempelvi)
endostatus_spike_temp_mods[[5]] <- glmmTMB(Spikelet ~ Endo*temp_cv + (1|Site/Plot) + (1|Plot/Population) , ziformula = ~1,family = nbinom2,data=demography_tempelvi)
AICtab(endostatus_spike_temp_mods)
summary(endostatus_spike_temp_mods[[4]])


coefspk<-fixef(endostatus_spike_temp_mods[[4]])
s.isEndo0<-coefspk[1]$cond[1]			## intercept for Endo 0
s.ssEndo0<-coefspk[1]$cond[3]			## slope for for Endo 0

s.isEndo1 <-coefspk[1]$cond[2]+ s.isEndo0	## intercept for Endo 1
s.ssEndo1 <- s.ssEndo0	## slope for Endo 1



temp_cv_seq<-seq(min(na.omit(demography_tempelvi$temp_cv)),max(na.omit(demography_tempelvi$temp_cv)),length.out=20)
temp_mean_seq<-seq(min(na.omit(demography_tempelvi$temp_mean)),max(na.omit(demography_tempelvi$temp_mean)),length.out=20)


pdf("/Users/jmoutouama/Dropbox/Miller Lab/github/ELVI-endo-density/Figure/Temperature.pdf",height=3,width =9,useDingbats = F)
# layout(mat = layout.matrix,
#        heights = rep(c(1, 1,1),2), # Heights of the two rows
#        widths = c(3, 3,3))
# layout.show(3)
# par(oma=c(5,1,0.5,0.5))

par(mfrow=c(1,3),mar=c(5,5,2,1))
with(demography_tempelvi,{
  # par(mar=c(4,4,2,4))
  plot(temp_cv,grow,cex.lab=1.5,xlab="Soil temperature (CV)",ylab="Growth rate");box()
  points(jitter(temp_cv),jitter(grow),bg=cbPalette,pch=21,col=NA,cex=2,lwd=2)
  title("A",adj=0)
  lines(temp_cv_seq,(s.igEndo0 + s.sgEndo0*temp_cv_seq),col=cbPalette[1],lwd=3,)
  lines(temp_cv_seq,(s.igEndo1 + s.sgEndo1*temp_cv_seq),col=cbPalette[2],lwd=3,)
})

with(demography_tempelvi,{
  # par(mar=c(4,4,2,4))
  plot(temp_cv,Inf_23,cex.lab=1.5,xlab="Soil temperature (CV)",ylab="# Inflorescences");box()
  points(jitter(temp_cv),jitter(Inf_23),bg=cbPalette,pch=21,col=NA,cex=2,lwd=2)
  title("B",adj=0)
  lines(temp_cv_seq,exp(s.iiEndo0 + s.siEndo0*temp_cv_seq),col=cbPalette[1],lwd=3,)
  lines(temp_cv_seq,exp(s.iiEndo1 + s.siEndo1*temp_cv_seq),col=cbPalette[2],lwd=3,)
})

with(demography_tempelvi,{
  # par(mar=c(4,4,2,4))
  plot(temp_cv,Spikelet,cex.lab=1.5,xlab="Soil temperature (CV)",ylab="# Spikelets");box()
  points(jitter(temp_cv),jitter(Spikelet),bg=cbPalette,pch=21,col=NA,cex=2,lwd=2)
  title("C",adj=0)
  lines(temp_cv_seq,exp(s.isEndo0 + s.ssEndo0*temp_cv_seq),col=cbPalette[1],lwd=3,)
  lines(temp_cv_seq,exp(s.isEndo1 + s.ssEndo1*temp_cv_seq),col=cbPalette[2],lwd=3,)
  legend("topright",legend=c("E-","E+"),lwd=2,col =cbPalette ,bty="n",cex=1)
})

dev.off()


endostatus_growth_soil <- list()
endostatus_growth_soil[[1]] <-lmer(grow ~ 1 + (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_soil[[2]] <- lmer(grow ~  Endo + (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_soil[[3]] <- lmer(grow ~ Endo + water_mean + (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_soil[[4]] <- lmer(grow ~ Endo*water_mean + (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_soil[[5]] <- lmer(grow ~ Endo + water_cv+ (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_soil[[6]] <- lmer(grow ~ Endo*water_cv+ (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_soil[[7]] <- lmer(grow ~  Endo + (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_soil[[8]] <- lmer(grow ~ Endo + temp_mean + (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_soil[[9]] <- lmer(grow ~ Endo*temp_mean + (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_soil[[10]] <- lmer(grow ~ Endo + temp_cv+ (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
endostatus_growth_soil[[11]] <- lmer(grow ~ Endo*temp_cv+ (1|Site/Plot) + (1|Population) , data=demography_tempelvi,na.action = na.omit)
AICtab(endostatus_growth_soil)
