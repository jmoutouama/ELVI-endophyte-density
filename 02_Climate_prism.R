# Project: 
# Purpose: Create variables that most accurately reflect the climate across study area.
# Note: Raster files are too large to provide in public repository. They are stored on a local machine
# Authors: Jacob Moutouama
# Date last modified (Y-M-D): 2024-08-03
rm(list = ls())
# load packages
library(tidyverse) # a suite of packages for wrangling and tidying data
library(prism)     # package to access and download climate data
library(raster)    # the climate data comes in raster files- this package helps process those
library(stringr)
library(magrittr)
library(readxl)# read excel data
library(ggsci) # package for color blind color in ggplot 2
library(corrplot)# visualize the correlation  
library(terra)
library(zoo)
library(SPEI)
library(smplot2)
# PRISM data ----
# First, set a file path where prism data will be stored
options(prism.path = '/Users/jm200/Documents/Prism Range limit/')
#get_prism_monthlys(type="ppt",years=1994:2024,mon=1:12,keepZip = TRUE)

# Grab the prism data and compile the files
climate_data <- prism_archive_ls() %>%  
  pd_stack(.)  
climate_crs<-climate_data@crs@projargs
# Convert these locations to format that can be matched to Prism climate data
read.csv("https://www.dropbox.com/scl/fi/1eu5lhkg5mx7roj3zd7g0/Study_site.csv?rlkey=tonb6sswc7zqf123ct06t64yp&dl=1", stringsAsFactors = F) %>% 
  arrange(latitude)->garden ## common garden populations
garden_sites <- as.data.frame(garden)
coordinates(garden_sites) <- c('longitude', 'latitude')
proj4string(garden_sites) <- CRS(climate_crs)
# Extract climatic data from the raster stack for those sites 
climate_garden <- data.frame(coordinates(garden_sites), 
                             garden_sites$site_code, 
                             extract(climate_data, garden_sites))
# Reshape data. Col 1:3 are lat, long, and site ID. Col 4:ncol are climate data
# Column headers include date and climate type info
climate_garden <- climate_garden %>% 
  gather(date, value, 4:ncol(climate_garden))
# The column header includes the date and data type, but also some other metadata that we don't need
# Here, I remove the extra info from the column header
climate_garden$date <- gsub('PRISM_', '', climate_garden$date) %>% 
  gsub('stable_4kmM3_', '', .) %>% 
  gsub('provisional_4kmM3_', '', .) %>%
  gsub('_bil', '', .)

# Split header into type (precipitation or temperature), year, and month
climate_garden <- separate(climate_garden, 'date', 
                           into = c('clim', 'YearMonth'), 
                           sep = '_')
climate_garden <- separate(climate_garden, 'YearMonth',
                           into = c('year', 'month'),
                           sep = 4)
# Reshape data-- make a separate column for temperature and precipitation
climate_garden <- unique(climate_garden)
climate_garden <- climate_garden %>% 
  spread(clim, value) %>%
  rename(lon = longitude, lat = latitude, site = garden_sites.site_code)

# Calculate the potential evapotranspiration (PET) based on thornthwaite method and compute the water balance (P - PET).
climate_garden %>% 
  filter(site=="SON") %>% 
  mutate(PET = thornthwaite(tmean,30.27347), BAL=ppt-PET)->climate_garden_SON
climate_garden %>% 
  filter(site=="KER") %>% 
  mutate(PET = thornthwaite(tmean,30.02844), BAL=ppt-PET)->climate_garden_KER
climate_garden %>% 
  filter(site=="BAS") %>% 
  mutate(PET = thornthwaite(tmean,30.09234), BAL=ppt-PET)->climate_garden_BAS
climate_garden %>% 
  filter(site=="BFL") %>% 
  mutate(PET = thornthwaite(tmean,30.28487), BAL=ppt-PET)->climate_garden_BFL
climate_garden %>% 
  filter(site=="LAF") %>% 
  mutate(PET = thornthwaite(tmean,30.30477), BAL=ppt-PET)->climate_garden_LAF
climate_garden %>% 
  filter(site=="COL") %>% 
  mutate(PET = thornthwaite(tmean,30.56930), BAL=ppt-PET)->climate_garden_COL
climate_garden %>% 
  filter(site=="HUN") %>% 
  mutate(PET = thornthwaite(tmean,30.74281), BAL=ppt-PET)->climate_garden_HUN

# Convert to a ts (time series) for convenience
climate_garden_SON_ts <- ts(climate_garden_SON[, -c(1,2,3,4,5)], start=c(1990,1), end = c(2024, 12), frequency = 12)
plot(climate_garden_SON_ts)
spei_SON <- SPEI::spei(climate_garden_SON_ts[, "BAL"], scale = 3)
plot(spei_SON)
spei_SON_values<-as.data.frame(spei_SON$fitted)
names(spei_SON_values)<-"SPEI"
climate_garden_SON_SPEI<-data.frame(climate_garden_SON,spei_SON_values)

climate_garden_KER_ts <- ts(climate_garden_KER[, -c(1,2,3,4,5)], start=c(1990, 01), end = c(2024, 12), frequency = 12)
plot(climate_garden_KER_ts)
spei_KER <- SPEI::spei(climate_garden_KER_ts[, "BAL"], scale = 3)
plot(spei_KER)
spei_KER_values<-as.data.frame(spei_KER$fitted)
names(spei_KER_values)<-"SPEI"
climate_garden_KER_SPEI<-data.frame(climate_garden_KER,spei_KER_values)

climate_garden_BAS_ts <- ts(climate_garden_BAS[, -c(1,2,3,4,5)], start=c(1990, 01), end = c(2024, 12), frequency = 12)
plot(climate_garden_BAS_ts)
spei_BAS <- SPEI::spei(climate_garden_BAS_ts[, "BAL"], scale = 3)
plot(spei_BAS)
spei_BAS_values<-as.data.frame(spei_BAS$fitted)
names(spei_BAS_values)<-"SPEI"
climate_garden_BAS_SPEI<-data.frame(climate_garden_BAS,spei_BAS_values)

climate_garden_BFL_ts <- ts(climate_garden_BFL[, -c(1,2,3,4,5)], start=c(1990, 01), end = c(2024, 12), frequency = 12)
plot(climate_garden_BFL_ts)
spei_BFL <- SPEI::spei(climate_garden_BFL_ts[, "BAL"], scale = 3)
plot(spei_BFL)
spei_BFL_values<-as.data.frame(spei_BFL$fitted)
names(spei_BFL_values)<-"SPEI"
climate_garden_BFL_SPEI<-data.frame(climate_garden_BFL,spei_BFL_values)

climate_garden_LAF_ts <- ts(climate_garden_LAF[, -c(1,2,3,4,5)], start=c(1990, 01), end = c(2024, 12), frequency = 12)
plot(climate_garden_LAF_ts)
spei_LAF <- SPEI::spei(climate_garden_LAF_ts[, "BAL"], scale = 3)
plot(spei_LAF)
spei_LAF_values<-as.data.frame(spei_LAF$fitted)
names(spei_LAF_values)<-"SPEI"
climate_garden_LAF_SPEI<-data.frame(climate_garden_LAF,spei_LAF_values)

climate_garden_COL_ts <- ts(climate_garden_COL[, -c(1,2,3,4,5)], start=c(1990, 01), end = c(2024, 12), frequency = 12)
plot(climate_garden_COL_ts)
spei_COL <- SPEI::spei(climate_garden_COL_ts[, "BAL"], scale = 3)
plot(spei_COL)
spei_COL_values<-as.data.frame(spei_COL$fitted)
names(spei_COL_values)<-"SPEI"
climate_garden_COL_SPEI<-data.frame(climate_garden_COL,spei_COL_values)

climate_garden_HUN_ts <- ts(climate_garden_HUN[, -c(1,2,3,4,5)], start=c(1990, 01), end = c(2024, 12), frequency = 12)
plot(climate_garden_HUN_ts)
spei_HUN <- SPEI::spei(climate_garden_HUN_ts[, "BAL"], scale = 3)
plot(spei_HUN)
spei_HUN_values<-as.data.frame(spei_HUN$fitted)
names(spei_HUN_values)<-"SPEI"
climate_garden_HUN_SPEI<-data.frame(climate_garden_HUN,spei_HUN_values)

# Put all the sites together
climate_garden_SPEI<-bind_rows(climate_garden_HUN_SPEI,
                               climate_garden_COL_SPEI,
                               climate_garden_LAF_SPEI,
                               climate_garden_BAS_SPEI,
                               climate_garden_BFL_SPEI,
                               climate_garden_KER_SPEI,
                               climate_garden_SON_SPEI)

# Subset the data collection period
climate_garden_SPEI %>% 
  unite("longdate",year:month,sep= "-") %>% 
  filter(longdate>as.yearmon("2023-05") & longdate<as.yearmon("2024-06"))->climate_garden_SPEI_2023_2024

## Plot the daily trend for temperature and soil moisture from start to end
climate_garden_SPEI_2023_2024 %>% 
  group_by(site) %>% 
  summarise(
    mean_ppt=mean(ppt),
    mean_temp=mean(tmean),
    mean_pet=mean(PET),
    mean_spei=mean(SPEI))->prism_means
prism_means<-as.data.frame(prism_means)
# saveRDS(prism_means, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Data/prism_means.rds')

#  Time scale provide insight into long-term drought trends.
# SPEI > 0: Wet conditions.
# SPEI < 0: Dry conditions (drought).
# SPEI ≤ -1.5: Moderate to severe drought.

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Climate_prism.pdf",width=14,height=10,useDingbats = F)
par(mar=c(5,5,2,3),mfrow=c(2,2))
barplot(prism_means[order(prism_means[,2],decreasing=FALSE),][,2],names.arg=prism_means[order(prism_means[,2],decreasing=FALSE),][,1],col="#E69F00",xlab="Sites", ylab="Mean",main="",ylim=c(0,200))
mtext("Precipitation",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "A",side = 3, adj = 0,cex=1.2)
barplot(prism_means[order(prism_means[,3],decreasing=FALSE),][,3],names.arg=prism_means[order(prism_means[,3],decreasing=FALSE),][,1],col="#E69F00",xlab="Sites", ylab="Mean",main="",ylim=c(0,25))
mtext("Temperature",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "B",side = 3, adj = 0,cex=1.2)
barplot(prism_means[order(prism_means[,4],decreasing=FALSE),][,4],names.arg=prism_means[order(prism_means[,4],decreasing=FALSE),][,1],col="#E69F00",xlab="Sites", ylab="Mean",main="")
mtext("Potential evapotranspiration (PET)",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "C",side = 3, adj = 0,cex=1.2)
barplot(prism_means[order(prism_means[,5],decreasing=FALSE),][,5],names.arg=prism_means[order(prism_means[,5],decreasing=FALSE),][,1],col="#E69F00",xlab="Sites", ylab="Mean",main="")
mtext("Standardized Precipitation Evapotranspiration Index",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "D",side = 3, adj = 0,cex=1.2)
dev.off()


site_names <- c("LAF"="Lafayette",
                "HUN"="Huntville",
                "KER"="Kerville",
                "BAS"="Bastrop",
                "COL"="College Station",
                "BFL" ="Brackenridge",
                "SON" ="Sonora")
climate_garden_SPEI_2023_2024 %>% 
  ggplot(aes(x=as.Date(as.yearmon(longdate)), y=ppt))+
  geom_line(aes(colour=site))+
  #ggtitle("d")+
  # scale_color_manual(values = cbp1)+
  # scale_fill_manual(values = cbp1)+
  theme_bw()+ 
  theme(legend.position = "none",
        axis.text.x = element_text(size=4.5,color="black", angle=0),
        plot.title =element_text(size=14, color="black",angle=0))+
  labs( y="Daily precipitation  (°C)", x="Month")+
  #ylim=c(0,100)+
  #facet_grid(~site,labeller = labeller(site=site_names))+
  facet_grid(~factor(site,levels=c("HUN","LAF","COL","BAS","BFL","KER","SON")))+
  geom_hline(data=prism_means,aes(yintercept = mean_ppt,colour=site))->figpptsite_prism

ggplot(data = climate_garden_SPEI_2023_2024, mapping = aes(x = tmean, y = ppt)) +
  geom_point(shape = 21, fill = "#0f993d", color = "white", size = 3) +
  sm_statCorr()
