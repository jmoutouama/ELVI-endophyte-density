# Project: 
# Purpose: Maps of Common garden experiment for Agrostis hyemalis, Elymus virginicus and Poa autumnalis across a climatic gradient. 
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
library(maptools)
# Climatic data----
## Data from PRISM---- 
# making a folder to store prism data
options(prism.path = "/Users/jm200/Documents/PRISM")
# getting monthly data for mean temp and precipitation
# takes a long time the first time, but can skip when you have raster files saved on your computer.
# get_prism_monthlys(type = "tmean", years = 1994:2024, mon = 1:12, keepZip = FALSE)
# get_prism_monthlys(type = "ppt", years = 1994:2024, mon = 1:12, keepZip = FALSE)
prism_set_dl_dir("/Users/jm200/Documents/PRISM")
prism_archive_ls()

# Common garden, natural  and source populations locations ---- 
read.csv("https://www.dropbox.com/scl/fi/1eu5lhkg5mx7roj3zd7g0/Study_site.csv?rlkey=tonb6sswc7zqf123ct06t64yp&dl=1", stringsAsFactors = F) %>% 
  dplyr::select(latitude,longitude) %>%
  unique() %>% 
  arrange(latitude)->garden ## common garden populations

read.csv("https://www.dropbox.com/scl/fi/go448bqe9z6meisgkd6g9/source_pop.csv?rlkey=oeutzzf4lgo616zeam06tul8t&dl=1", stringsAsFactors = F) %>% 
  dplyr::select(latitude,longitude) %>% 
  unique() %>% 
  arrange(latitude)->source ## source populations

# Agrostis hyemalis-----
#aghy_occ_raw <- gbif(genus="Agrostis",species="hyemalis",download=TRUE) 
#saveRDS(aghy_occ_raw, file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/ELVI Model output/occurence/aghy_occ_raw.rds")
aghy_occ_raw<-readRDS(url("https://www.dropbox.com/scl/fi/ijl1i7964qxzcxvm7blz8/aghy_occ_raw.rds?rlkey=msbs3xzkjc719uld8yw7hb6cj&dl=1"))
head(aghy_occ_raw) 

aghy_occ_raw %>% 
  filter(!is.na(lat) & !is.na(lon) & country=="United States") %>% 
  unique() %>% 
  dplyr::select(country,lon, lat,year)%>% 
  arrange(lat)->aghy

# Elymus virginicus----
#elvi_occ_raw <- gbif(genus="Elymus",species="virginicus",download=TRUE) 
#saveRDS(elvi_occ_raw, file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/ELVI Model output/occurence/elvi_occ_raw.rds")
elvi_occ_raw<-readRDS(url("https://www.dropbox.com/scl/fi/0ssa5gepxyz28b7ykw1x8/elvi_occ_raw.rds?rlkey=4dx0q4lw2112droh73hmh7xte&dl=1"))
head(elvi_occ_raw) 
elvi_occ_raw %>% 
  filter(!is.na(lat) & !is.na(lon) &  as.numeric(coordinateUncertaintyInMeters <=100) & country=="United States") %>% 
  unique() %>% 
  dplyr::select(country,lon, lat,year)%>% 
  arrange(lat)->elvi

# Poa autumnalis ---
#poa_occ_raw <- gbif(genus="Poa",species="autumnalis",download=TRUE) 
#saveRDS(poa_occ_raw, file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/ELVI Model output/occurence/poa_occ_raw.rds")
poa_occ_raw<-readRDS(url("https://www.dropbox.com/scl/fi/oip7ndyf0d99rqxcqxb0q/poa_occ_raw.rds?rlkey=920uql1gd4gahnh8utw9fz96l&dl=1"))
head(poa_occ_raw) 
poa_occ_raw %>% 
  filter(!is.na(lat) & !is.na(lon) & as.numeric(coordinateUncertaintyInMeters <=100)  & country=="United States") %>% 
  unique() %>% 
  dplyr::select(country,lon, lat,year)%>% 
  arrange(lat)->poa

# Georeferencing the occurences -----
coordinates(aghy) <- ~ lon + lat
coordinates(elvi) <- ~ lon + lat
coordinates(poa) <- ~ lon + lat
coordinates(garden) <- ~ longitude + latitude
coordinates(source) <- ~ longitude + latitude
CRS1 <- CRS("+init=epsg:4326") # WGS 84
crs(aghy) <- CRS1
crs(poa) <- CRS1
crs(elvi) <- CRS1
crs(garden) <- CRS1
crs(source) <- CRS1

# Study area shapefile ----
US<-terra::vect("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/data/USA_vector_polygon/States_shapefile.shp")
plot(US)
aghy_d <- US[(US$State_Name %in% c("TEXAS","LOUISIANA","ARKANSAS","OKLAHOMA","MISSISSIPPI","ALABAMA","TENNESSEE","KENTUCKY","GEORGIA","NORTH CAROLINA","SOUTH CAROLINA","MISSOURI","KANSAS","WEST VIRGINIA","VIRGINIA","ILLINOIS","PENNSYLVANIA","INDIANA","OHIO","MARYLAND","NEW JERSEY"  )), ]
elvi_d <- US[(US$State_Name %in% c("TEXAS","LOUISIANA","ARKANSAS","OKLAHOMA","MISSISSIPPI","ALABAMA","TENNESSEE","KENTUCKY","GEORGIA","NORTH CAROLINA","SOUTH CAROLINA","MISSOURI","KANSAS","WEST VIRGINIA","VIRGINIA","ILLINOIS","PENNSYLVANIA","INDIANA","OHIO","MARYLAND","NEBRASKA","NEW JERSEY","IOWA","WISCONSIN","MICHIGAN","NEW YORK","SOUTH DAKOTA","RHODE ISLAND","VERMONT","NEW HAMPSHIRE","MASSACHUSETTS","NORTH DAKOTA","MINNESOTA","CONNECTICUT","NEW MEXICO","COLORADO","FLORIDA","DELAWARE","MARYLAND")), ]
poa_d <- US[(US$State_Name %in% c("TEXAS","LOUISIANA","ARKANSAS","OKLAHOMA","MISSISSIPPI","ALABAMA","TENNESSEE","KENTUCKY","GEORGIA","NORTH CAROLINA","SOUTH CAROLINA","MISSOURI","VIRGINIA","MARYLAND")), ]

plot(aghy_d)
plot(elvi_d)
# calculating the cumulative precipitation for each year and for each season within the year
ppt_annual <- list()
for(y in 1994:2024){
  ppt_annual[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y))))
}
# Taking the mean of the cumulative precipitation values
ppt_annual_norm <- terra::mean(terra::rast(unlist(ppt_annual)))
crs(ppt_annual_norm)<-CRS1
crop_ppt_annual_aghy <- terra::crop(ppt_annual_norm, US,mask=TRUE)
crop_ppt_annual_elvi <- terra::crop(ppt_annual_norm, elvi_d,mask=TRUE)
crop_ppt_annual_poa <- terra::crop(ppt_annual_norm, poa_d,mask=TRUE)

col_precip <- terrain.colors(30)

# Maps (Figure 1) ----
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/clim_map.pdf",width=9,height=8)
op<-par(mfrow = c(2,2), mar=c(0,1,3.75,1),oma = c(0, 1, 1, 0))
#par(mar=c(5,0,4,0),mfrow=c(2,2))
plot(crop_ppt_annual_aghy,xlab="Longitude",ylab="Latitude",col=col_precip,cex.lab=1.2)
#text(x=-85, y=34, "Precipitation (mm)", srt=-90, cex=0.8, xpd=NA, pos=4)
plot(aghy_d,add=T)
plot(aghy,add=T,pch = 23,col="grey50",bg="grey",cex =0.55)
plot(garden,add=T,pch = 3,col="black",cex =2)
plot(source,add=T,pch = 21,col="black",bg="red",cex =1)
mtext(~ italic("Agrostis hyemalis"),side = 3, adj = 0.5,cex=1.25,line=0.2)

plot(crop_ppt_annual_elvi,xlab="Longitude",ylab="",col=col_precip,cex.lab=1.2)
#text(x=-85, y=34, "Precipitation (mm)", srt=-90, cex=0.8, xpd=NA, pos=4)
plot(elvi_d,add=T)
plot(elvi,add=T,pch = 23,col="grey50",bg="grey",cex =0.55)
plot(garden,add=T,pch = 3,col="black",cex =2)
plot(source,add=T,pch = 21,col="black",bg="red",cex =1)
mtext(~ italic ("Elymus virginicus"),side = 3, adj = 0.5,cex=1.25,line=0.2)

par(mar=c(0,3,3.75,1))
plot(crop_ppt_annual_poa,xlab="Longitude",ylab="Latitude",col=col_precip,cex.lab=1.2)
#text(x=-85, y=34, "Precipitation (mm)", srt=-90, cex=0.8, xpd=NA, pos=4)
plot(poa_d,add=T)
plot(poa,add=T,pch = 23,col="grey50",bg="grey",cex =0.55)
plot(garden,add=T,pch = 3,col="black",cex =2)
plot(source,add=T,pch = 21,col="black",bg="red",cex =1)
mtext( ~ italic ("Poa autumnalis"),side = 3, adj = 0.5,cex=1.25,line=0.2)
legend(-106, 28.5, 
       legend=c( "GBIF occurences","Common garden sites","Source populations"),
       pch = c(23,3,21),
       pt.cex=c(0.55,1,1),
       col = c("grey50","black","black"),
       pt.bg=c("grey","black","red"),
       cex = 0.7, 
       bty = "n", 
       horiz = F , 
)
par(op)
dev.off()





