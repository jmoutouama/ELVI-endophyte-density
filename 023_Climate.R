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
# Extract the  extracted from the raster stack for those sites 
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

# Make year, month and day numeric variables
climate_garden$year  <- as.numeric(climate_garden$year)
climate_garden$month  <- as.numeric(climate_garden$month)
# Order data by site
climate_garden <- climate_garden[order(climate_garden$site),]
climate_garden<-climate_garden %>% 
  unite("longdate",year:month,sep= "-")

## to standardize census data with prism data
climate_garden %>% 
  filter(site==c("BAS","BFL","COL","HUN","LAF"))->climate_garden_sites

## Plot the daily trend for temperature and soil moisture from start to end
climate_garden_sites %>% 
  group_by(site) %>% 
  summarise(
            mean_ppt=mean(ppt))->prism_means
prism_means<-as.data.frame(prism_means)
barplot(prism_means[order(prism_means[,2],decreasing=FALSE),][,2],names.arg=prism_means[order(prism_means[,1],decreasing=FALSE),][,1],col="#E69F00",xlab="Precipitation", ylab="Mean",main="")


climate_garden_sites %>% 
  ggplot(aes(x=as.Date(as.yearmon(longdate)), y=ppt))+
  geom_line(aes(colour=site))+
  ggtitle("d")+
  # scale_color_manual(values = cbp1)+
  # scale_fill_manual(values = cbp1)+
  theme_bw()+ 
  theme(legend.position = "none",
        axis.text.x = element_text(size=4.5,color="black", angle=0),
        plot.title =element_text(size=14, color="black",angle=0))+
  labs( y="Daily precipitation  (°C)", x="Month")+
  #ylim=c(0,100)+
  #facet_grid(~Site,labeller = labeller(Site=site_names))+
  facet_grid(~factor(site,levels=c("LAF","HUN","COL","BAS","BFL")))+
  geom_hline(data=prism_means,aes(yintercept = mean_ppt,colour=site))->figpptsite_prism



