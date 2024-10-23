# Project: 
# Purpose: Maps of Common garden experiment for Poa arachnifera across a climatic gradient. 
# Authors: Jacob Moutouama
# Date last modified (Y-M-D): 

# remove all objects and clear workspace
rm(list = ls(all=TRUE))

# load packages
library(sp)
library(rgdal)
library(raster)
library(terra)
library(maptools)
library(rgeos)
library(RColorBrewer)
library(tidyverse)
library(dismo)
library(prism)
# Climatic data----
## Data from PRISM---- 
# making a folder to store prism data
options(prism.path = "/Users/jm200/Documents/PRISM")
# getting monthly data for mean temp and precipitation
# takes a long time the first time, but can skip when you have raster files saved on your computer.
get_prism_monthlys(type = "tmean", years = 1990:2020, mon = 1:12, keepZip = FALSE)
get_prism_monthlys(type = "ppt", years = 1990:2020, mon = 1:12, keepZip = FALSE)

prism_set_dl_dir("/Users/jm200/Documents/PRISM")
prism_archive_ls()
# Calculating mean, standard deviation and CV in temp for each year




temp_data <- ppt_data<- list()

for(y in 1990:2020){
  temp_data <- prism_stack(prism_archive_ls()())
  print(y)
}





for(y in 1990:2020){
  temp_annual_mean[[y]] <- mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = y))))
  temp_annual_sd[[y]] <- stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = y))))
  print(y)
  #temp_annual_cv[[y]] <-   temp_annual_mean[[y]]/temp_annual_cv[[y]]
}

plot(temp_annual_mean)

plot(tmean_annual_norm)



# Common garden, natural  and source populations locations ---- 
read.csv("https://www.dropbox.com/scl/fi/1eu5lhkg5mx7roj3zd7g0/Study_site.csv?rlkey=tonb6sswc7zqf123ct06t64yp&dl=1", stringsAsFactors = F) %>% 
  dplyr::select(latitude,longitude) %>%
  unique() %>% 
  arrange(latitude)->garden ## common garden populations

read.csv("https://www.dropbox.com/scl/fi/go448bqe9z6meisgkd6g9/source_pop.csv?rlkey=oeutzzf4lgo616zeam06tul8t&dl=1", stringsAsFactors = F) %>% 
  dplyr::select(latitude,longitude) %>% 
  unique() %>% 
  arrange(latitude)->source ## source populations

#dir.create("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Data/occurence")
elvi_occ_raw <- gbif(genus="Elymus",species="virginicus",download=TRUE) 
head(elvi_occ_raw) 
elvi_occ <- subset(elvi_occ_raw,(!is.na(lat))&(!is.na(lon))) # here we remove erroneous coordinates, where either the latitude or longitude is missing
names(elvi_occ)
elvi_occ %>% 
  dplyr::select(country,lon, lat,year)%>% 
  dplyr::rename(longitude=lon,latitude=lat) %>% 
  filter(year %in% (1901:2024) & as.numeric(longitude >=-126.374160) &  as.numeric(longitude <=-95) & country=="United States") %>% 
  unique() %>% 
  arrange(latitude)->gbif

coordinates(gbif) <- ~ longitude + latitude
CRS1 <- CRS("+init=epsg:4326") # WGS 84
crs(gbif) <- CRS1
coordinates(garden) <- ~ longitude + latitude
crs(garden) <- CRS1
coordinates(source) <- ~ longitude + latitude
crs(source) <- CRS1

# Study area shapefile ----
study_area<-terra::vect("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/data/USA_vector_polygon/States_shapefile.shp")
study_area <- study_area[(study_area$State_Name %in% c("TEXAS","LOUISIANA")), ]
plot(study_area)

# Clip the climatic rasters
pptdorm_raster<-terra::rast(clim_1990_2019[[1]])
tempdorm_raster<-terra::rast(clim_1990_2019[[2]])
pptgrow_raster<-terra::rast(clim_1990_2019[[3]])
tempgrow_raster<-terra::rast(clim_1990_2019[[4]])
crop_pptdorm_raster <- terra::crop(pptdorm_raster, study_area,mask=TRUE)
crop_tempdorm_raster <- terra::crop(tempdorm_raster, study_area,mask=TRUE)
crop_pptgrow_raster <- terra::crop(pptgrow_raster, study_area,mask=TRUE)
crop_tempgrow_raster <- terra::crop(tempgrow_raster, study_area,mask=TRUE)


# Maps (Figure 1) ----
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/clim_map_v1.pdf",width=11,height=5)

par(mfrow=c(1,3))
plot(study_area,xlab="Longitude",ylab="Latitude",cex.lab=1.2)
plot(gbif,add=T,pch = 23,col="grey50",bg="grey",cex =0.55)
plot(garden,add=T,pch = 3,col="black",cex =2,bg=unique(poar_2015_2016$site_col))
plot(source,add=T,pch = 21,col="black",bg="red",cex =1)
mtext( "A",side = 3, adj = 0,cex=1.25,line=0.2)
legend(-106, 28.5, 
       legend=c( "GBIF occurences","Common garden sites","Source populations"),
       pch = c(23,3,21),
       pt.cex=c(0.55,2,1),
       col = c("grey50","black","black"),
       pt.bg=c("grey","black","red"),
       cex = 0.9, 
       bty = "n", 
       horiz = F , 
)

dev.off()







