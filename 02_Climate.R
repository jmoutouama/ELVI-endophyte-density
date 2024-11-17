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
# Path to acess the data
jacob_path <- "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density"
#jacob_path<-"D:/ELVI-endophyte-density"
# tom_path<-"C:/Users/tm9/Dropbox/github/ELVI-endophyte-density"
choose_path <- jacob_path
# Demographic data -----
# Merge the demographic census
datini<-read_csv(paste0(choose_path,"/Data/Initialdata.csv"))
dat23<-read_csv(paste0(choose_path,"/Data/census2023.csv"))
dat24<-read_csv(paste0(choose_path,"/Data/census2024.csv"))
datherbivory<-read_csv(paste0(choose_path,"/Data/herbivory.csv"))
# calculate the total spikelet for each census
dat23 %>%
  mutate(
    spi_mean_23 = rowMeans(across(Spikelet_A:Spikelet_C), na.rm = TRUE),
         spikelet_23 = ceiling(spi_mean_23)
    ) -> dat23_spike

dat24 %>%
  mutate(
    spi_mean_24 = rowMeans(across(Spikelet_A:Spikelet_C), na.rm = TRUE),
    spikelet_24 = ceiling(spi_mean_24),
    Inf_sum_24 = rowSums(across(attachedInf_24:brokenInf_24), na.rm =TRUE),
    Inf_24 = ceiling(Inf_sum_24)
  ) -> dat24_spike

## Merge the initial data with the 23 data and the 23 data with the 24 -----
datini23 <- left_join(x = datini,y =dat23_spike,by=c("Tag_ID"))
#names(datini23)
dat2324 <- left_join(x = datini23 ,y =dat24_spike,by=c("Tag_ID")) 
#names(dat2324)
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
#unique(dat2324_t_t1_herb)
#head(dat2324_t_t1_herb)
#view(dat2324_t_t1_herb)
# ELVLI data
dat2324_t_t1_herb %>% 
  filter(Species=="ELVI")->dat2324_t_t1_herb_ELVI
#view(dat2324_t_t1_herb_ELVI)
## Consider only 5 sites 
dat2324_t_t1_herb_ELVI %>% 
  filter(Site %in% c("BAS","BFL","COL","HUN","LAF","SON"))->dat2324_t_t1_herb_ELVI_clean
#view(dat2324_t_t1_herb_ELVI_clean)
#dim(dat2324_t_t1_herb_ELVI_clean)

## Find the starting and ending dates are correct
dat2324_t_t1_herb_ELVI_clean %>% 
  dplyr::select(Site,Species,date_23,date_24) %>% 
  group_by(Site) %>% 
  unique() %>% 
  mutate(duration=difftime(mdy(date_24),mdy(date_23),units = "days"))->census_dates

# HOBO data ----
## format date and separate year-month-day
list.files(
  path = paste0(choose_path, "/Data/HOBO data/"),
  pattern = "*.xlsx",
  full.names = TRUE
) %>% # Identify all excel files
  lapply(read_excel) %>%                              # Store all files in list
  bind_rows -> hobo_data_raw # get HOBO data

tidyr::separate(hobo_data_raw,
                "date",
                into = c('longdate', 'time'),
                sep = ' ') %>%
  tidyr::separate(
    'longdate',
    # Separate the ‘longdate’ column into separate columns for month, day and year using the separate() function.
    into = c('year', 'month', 'day'),
    sep = '-',
    remove = FALSE
  ) -> hobo_data_full 
#names(hobo_data_full)
## double check if the starting and ending dates are correct
hobo_data_full %>% 
  group_by(site) %>% 
  summarise(start=range(longdate)[1],
            end=range(longdate)[2],
            duration=as.Date(end)-as.Date(start))->hobo_dates 

## average over days to look at overall trend across sites
hobo_data_full %>%
  group_by(longdate, site, day) %>%
  summarise(daily_mean_moist = mean(water),
            daily_mean_temp = mean(temperature)) -> HOBO_daily

##Standardize census data with hobo data
HOBO_daily %>% 
  filter(site=="BAS" & longdate>as.Date("2023-06-22") & longdate<as.Date("2024-06-14"))->HOBO_BAS
HOBO_daily %>% 
  filter(site=="BFL" & longdate>as.Date("2023-06-23") & longdate<as.Date("2024-06-14"))->HOBO_BFL
HOBO_daily %>% 
  filter(site=="COL" & longdate>as.Date("2023-06-09") & longdate<as.Date("2024-06-25"))->HOBO_COL
HOBO_daily %>% 
  filter(site=="HUN" & longdate>as.Date("2023-06-07") & longdate<as.Date("2024-06-04"))->HOBO_HUN
HOBO_daily %>% 
  filter(site=="LAF" & longdate>as.Date("2023-06-13") & longdate<as.Date("2024-06-06"))->HOBO_LAF
HOBO_daily %>% 
  filter(site=="SON" & longdate>as.Date("2023-06-28") & longdate<as.Date("2024-06-28"))->HOBO_SON

HOBO_daily_all_sites<-rbind(HOBO_BAS,HOBO_BFL,HOBO_COL,HOBO_HUN,HOBO_LAF,HOBO_SON)
#unique(HOBO_daily_all_sites$site)
## Plot the daily trend for temperature and soil moisture from start to end
HOBO_daily_all_sites %>% 
  group_by(site) %>% 
  summarise(mean_temp=mean(daily_mean_temp),
            mean_moisture=mean(daily_mean_moist))->hobo_means

data_plotclim<-data.frame(site=c(HOBO_daily_all_sites$site,HOBO_daily_all_sites$site),daily_mean_clim=c(HOBO_daily_all_sites$daily_mean_temp,HOBO_daily_all_sites$daily_mean_moist),date=c(HOBO_daily_all_sites$longdate,HOBO_daily_all_sites$longdate),clim=c(rep("temp",nrow(HOBO_daily_all_sites)),rep("water",nrow(HOBO_daily_all_sites))))
#HOBO_daily_all_sites$site <- factor(HOBO_daily_all_sites$site,levels = c("LAF", "HUN", "BAS", "COL", "BFL"))

site_names <- c("LAF"="Lafayette",
                 "HUN"="Huntville",
                 "BAS"="Bastrop",
                "COL"="College Station",
                "BFL" ="Brackenridge")

cbp1 <- c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7")



HOBO_daily_all_sites %>% 
  #mutate(site = ordered(site, levels=c("BFL","BAS","COL","HUN","LAF"))) %>%
  ggplot(aes(x=as.Date(longdate, format= "%Y - %m - %d"), y=daily_mean_temp))+
  geom_line(aes(colour=site))+
  ggtitle("a")+
  # scale_color_manual(values = cbp1)+
  # scale_fill_manual(values = cbp1)+
  theme_bw()+ 
  theme(legend.position = "none",
        axis.text.x = element_text(size=4.5,color="black", angle=0),
        plot.title =element_text(size=14, color="black",angle=0))+
  labs( y="Daily soil temperature  (°C)", x="")+
  facet_grid(~site,labeller = labeller(site=site_names))+
  geom_hline(data=hobo_means,aes(yintercept = mean_temp,colour=site))->figtempsite



HOBO_daily_all_sites %>%   
  ggplot(aes(x=as.Date(longdate, format= "%Y - %m - %d"), y=daily_mean_moist))+
  geom_line(aes(colour=site))+
  ggtitle("c")+
  # scale_color_manual(values = cbp1)+
  # scale_fill_manual(values = cbp1)+
  theme_bw()+ 
  theme(legend.position = "none",
        axis.text.x = element_text(size=4.55, color="black",angle=0),
        plot.title =element_text(size=14, color="black",angle=0))+
  labs( y="Daily soil moisture (wfv)", x="")+
  facet_grid(~site,labeller = labeller(site=site_names))+
  geom_hline(data=hobo_means,aes(yintercept = mean_moisture,colour=site))->figmoistsite

# Prism data ----
# First, set a file path where prism data will be stored
options(prism.path = '/Users/jm200/Documents/Prism ELVI/')
# get_prism_dailys(type = "tmean",minDate="2023-01-10",maxDate = "2024-07-31", keepZip = TRUE)
# get_prism_dailys(type = "ppt",minDate="2023-01-10",maxDate = "2024-07-31", keepZip = TRUE)

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
  gsub('stable_4kmD2_', '', .) %>% 
  gsub('provisional_4kmD2_', '', .) %>%
  gsub('_bil', '', .)

# Split header into type (precipitation or temperature), year, and month
climate_garden <- separate(climate_garden, 'date', 
                           into = c('clim', 'YearMonthDay'), 
                           sep = '_')
climate_garden <- separate(climate_garden, 'YearMonthDay',
                           into = c('year', 'MonthDay'),
                           sep = 4)
climate_garden <- separate(climate_garden, 'MonthDay',
                           into = c('month', 'day'),
                           sep = 2)
# Reshape data-- make a separate column for temperature and precipitation
climate_garden <- unique(climate_garden)
climate_garden <- climate_garden %>% 
  spread(clim, value) %>%
  rename(lon = longitude, lat = latitude, Site = garden_sites.site_code)

# Make year, month and day numeric variables
climate_garden$year  <- as.numeric(climate_garden$year)
climate_garden$month  <- as.numeric(climate_garden$month)
climate_garden$day  <- as.numeric(climate_garden$day)
# Order data by LTER site
climate_garden <- climate_garden[order(climate_garden$Site),]
climate_garden<-climate_garden %>% 
  unite("longdate",year:day,sep= "-")

## to standardize census data with prism data
climate_garden %>% 
  filter(Site=="BAS" & longdate>as.Date("2023-06-22") & longdate<as.Date("2024-06-14"))->climate_garden_BAS
climate_garden %>% 
  filter(Site=="BFL" & longdate>as.Date("2023-06-23") & longdate<as.Date("2024-06-14"))->climate_garden_BFL
climate_garden %>% 
  filter(Site=="COL" & longdate>as.Date("2023-06-09") & longdate<as.Date("2024-06-25"))->climate_garden_COL
climate_garden %>% 
  filter(Site=="HUN" & longdate>as.Date("2023-06-07") & longdate<as.Date("2024-06-04"))->climate_garden_HUN
climate_garden %>% 
  filter(Site=="LAF" & longdate>as.Date("2023-06-13") & longdate<as.Date("2024-06-06"))->climate_garden_LAF

climate_garden_daily_prism<-rbind(climate_garden_BAS,climate_garden_BFL,climate_garden_COL,climate_garden_HUN,climate_garden_LAF)

## Plot the daily trend for temperature and soil moisture from start to end
climate_garden_daily_prism %>% 
  group_by(Site) %>% 
  summarise(mean_temp=mean(tmean),
            mean_ppt=mean(ppt))->prism_means

climate_garden_daily_prism %>% 
  #mutate(site = ordered(site, levels=c("BFL","BAS","COL","HUN","LAF"))) %>%
  ggplot(aes(x=as.Date(longdate, format= "%Y - %m - %d"), y=tmean))+
  geom_line(aes(colour=Site))+
  ggtitle("b")+
  # scale_color_manual(values = cbp1)+
  # scale_fill_manual(values = cbp1)+
  theme_bw()+ 
  theme(legend.position = "none",
        axis.text.x = element_text(size=4.5,color="black", angle=0),
        plot.title =element_text(size=14, color="black",angle=0))+
  labs( y="Daily air temperature  (°C)", x="")+
  facet_grid(~Site,labeller = labeller(Site=site_names))+
  geom_hline(data=prism_means,aes(yintercept = mean_temp,colour=Site))->figtempsite_prism

climate_garden_daily_prism %>% 
  #mutate(site = ordered(site, levels=c("BFL","BAS","COL","HUN","LAF"))) %>%
  ggplot(aes(x=as.Date(longdate, format= "%Y - %m - %d"), y=ppt))+
  geom_line(aes(colour=Site))+
  ggtitle("d")+
  # scale_color_manual(values = cbp1)+
  # scale_fill_manual(values = cbp1)+
  theme_bw()+ 
  theme(legend.position = "none",
        axis.text.x = element_text(size=4.5,color="black", angle=0),
        plot.title =element_text(size=14, color="black",angle=0))+
  labs( y="Daily precipitation  (°C)", x="Month")+
  facet_grid(~Site,labeller = labeller(Site=site_names))+
  geom_hline(data=prism_means,aes(yintercept = mean_temp,colour=Site))->figpptsite_prism

# Regression to find the relationship between prism data and HOBO data
dat_reg<-data.frame(climate_garden_daily_prism,HOBO_daily_all_sites)
dat_reg<-dat_reg[,-c(1,2,3,4,7,8,9)]
corr_dat_reg <- cor(dat_reg)
library(corrplot)
corrplot(corr_dat_reg, tl.col = "brown", tl.srt = 45, bg = "White",
         title = "",
         type = "lower")

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Figregresion_temp.pdf",height =5,width=5,useDingbats = F)
(dat_reg %>% 
  ggplot(aes(tmean,daily_mean_temp))+
  geom_point()+
  geom_smooth(method = "lm",se=TRUE,fill = "blue", color = "blue")+
  labs( y="Soil temperature (°C)", x="Air temperature (°C)")+
  theme_bw()+
  theme(axis.text.x = element_text(size=4.5,color="black", angle=0))->Figregresion_temp)
dev.off()

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/climatesite.pdf",height =10,width=9,useDingbats = F)
(Figclimatesite<-ggpubr::ggarrange(figtempsite,figmoistsite,figtempsite_prism,figpptsite_prism,common.legend = FALSE,ncol = 1, nrow = 4))
dev.off()

