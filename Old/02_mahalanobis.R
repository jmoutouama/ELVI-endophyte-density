# Ecological niche modelling and niche centrality of Epichloe endophyte host species 
rm(list = ls())
# Load packages
#devtools::install_github('luismurao/ntbox')
library("prism")
library("raster")
library("dismo")
library("usdm")
library("car")
library("FactoMineR")
library("factoextra")
library("corrplot")
library("ENMTools")
library("vip")
library("pdp")
library("fastshap")
library("CalibratR")
library("maptools")
library("rgeos")
library("leaflet")
library("tidyverse")
# install.extras()
library("geodata")
library("terra")
library("rgdal")
library("sp")
library("ENMeval")
library("mecofun")
library("mgcv")
library("randomForest")
library("spThin")
library("ntbox")
# if( !("rJava" %in% rownames(installed.packages()))  ){
#   install.packages("rJava",repos="http://cran.r-project.org")
# }

# if(Sys.info()["sysname"] != "Windows" ){
#   dyn.load('/Library/Java/JavaVirtualMachines/temurin-17.jdk/Contents/Home/lib/server/libjvm.dylib')
# }
library("rJava")
library(geosphere)
library(oce)
library(sf)
library(ggplot2)
set.seed(13)
# Study area shapefile ----
study_area<-terra::vect("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/data/USA_vector_polygon/States_shapefile.shp")
study_area <- study_area[(study_area$State_Name %in% c("TEXAS","LOUISIANA")), ]
# Climatic data----
## Data from PRISM---- 
# making a folder to store prism data
options(prism.path = '/Users/jm200/Documents/Prism Range limit/')
# getting monthly data for mean temp and precipitation
# takes a long time the first time, but can skip when you have raster files saved on your computer.
# get_prism_monthlys(type = "tmean", years = 1994:2024, mon = 1:12, keepZip = FALSE)
# get_prism_monthlys(type = "ppt", years = 1994:2024, mon = 1:12, keepZip = FALSE)

# pulling out values to get normals for old and new time periods
tmean_annual_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1994:2024))))
tmean_spring_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1994:2024, mon = 1:4))))
tmean_summer_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1994:2024, mon = 5:8))))
tmean_autumn_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1994:2024, mon = 9:12))))

# calculating standard deviation in temp
tmean_annual_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1994:2024))))
tmean_spring_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1994:2024, mon = 1:4))))
tmean_summer_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1994:2024, mon = 5:8))))
tmean_autumn_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1994:2024, mon = 9:12))))

# calculating the cumulative precipitation for each year and for each season within the year
ppt_annual <- ppt_spring <- ppt_summer <- ppt_autumn <- ppt_winter<- list()
for(y in 1994:2024){
  ppt_annual[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y))))
  ppt_spring[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y, mon = 1:4))))
  ppt_summer[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y, mon = 5:8))))
  ppt_autumn[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y, mon = 9:12)))) 
}

# Taking the mean of the cumulative precipation values
ppt_annual_norm <- terra::mean(terra::rast(unlist(ppt_annual)))
ppt_spring_norm <- terra::mean(terra::rast(unlist(ppt_spring)))
ppt_summer_norm <- terra::mean(terra::rast(unlist(ppt_summer)))
ppt_autumn_norm <- terra::mean(terra::rast(unlist(ppt_autumn)))

#calculating the standard devation in precip
ppt_annual_sd <- terra::stdev(terra::rast(unlist(ppt_annual)))
ppt_spring_sd <- terra::stdev(terra::rast(unlist(ppt_spring)))
ppt_summer_sd <- terra::stdev(terra::rast(unlist(ppt_summer)))
ppt_autumn_sd <- terra::stdev(terra::rast(unlist(ppt_autumn)))

## Variance inflation factor (VIF)  to have  a measure of multicollinearity among the  variables----  
US_worldclim_norm<-terra::rast(list(tmean_spring_norm,tmean_summer_norm,tmean_autumn_norm,
                                           tmean_spring_sd,tmean_summer_sd,tmean_autumn_sd,
                                           ppt_spring_norm,ppt_summer_norm,ppt_autumn_norm,
                                           ppt_spring_sd,ppt_summer_sd,ppt_autumn_sd))
names(US_worldclim_norm) <- c("tmean_spring_norm","tmean_summer_norm","tmean_autumn_norm",
                              "tmean_spring_sd","tmean_summer_sd","tmean_autumn_sd",
                              "ppt_spring_norm","ppt_summer_norm","ppt_autumn_norm",
                              "ppt_spring_sd","ppt_summer_sd","ppt_autumn_sd")
US_worldclim_norm_stack<-stack(US_worldclim_norm)
plot(US_worldclim_norm_stack)

# Agrostis hyemalis-----
aghy_occ_raw <- gbif(genus="Agrostis",species="hyemalis",download=TRUE) 
head(aghy_occ_raw) 
#unique(aghy_occ_raw$basisOfRecord)
# aghy_occ_raw %>% 
#   filter(!is.na(lat) & !is.na(lon) & basisOfRecord=="HUMAN_OBSERVATION")->aghy_occ

aghy_occ <- subset(aghy_occ_raw,(!is.na(lat))&(!is.na(lon))&basisOfRecord==HUMAN_OBSERVATION) # here we remove erroneous coordinates, where either the latitude or longitude is missing
aghy_occ %>% 
  dplyr::select(country,lon, lat,year)%>% 
  dplyr::rename(longitude=lon,latitude=lat) %>% 
  filter(year %in% (1901:2024) & country=="United States") %>% 
  unique() %>% 
  arrange(latitude)->aghy

aghy_occ %>% 
  dplyr::select(country,lon, lat,year)%>% 
  dplyr::rename(longitude=lon,latitude=lat) %>% 
  filter(year %in% (1901:2024) & as.numeric(longitude >=-94) &  as.numeric(longitude <=-92) & as.numeric(latitude >=29.5) &  as.numeric(latitude <=32.5) & country=="United States") %>% 
  unique() %>% 
  arrange(latitude)->aghy2

aghy_occ %>% 
  dplyr::select(country,lon, lat,year)%>% 
  dplyr::rename(longitude=lon,latitude=lat) %>% 
  filter(year %in% (1901:2024) & as.numeric(longitude >=-92) &  as.numeric(longitude <=-89.5) & as.numeric(latitude >=29.5) &  as.numeric(latitude <=31) & country=="United States") %>% 
  unique() %>% 
  arrange(latitude)->aghy3

aghy<-rbind(aghy1,aghy2,aghy3)
class(aghy)
# Centroid
# Plot the polygon
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
centroid_aghy<-centroid(aghy[,2:3])
aghy_hull <- aghy[chull(aghy[, c('longitude', 'latitude')]),]
aghy_hull <- rbind(aghy_hull, aghy_hull[1,])
plot(study_area)
polygon(aghy_hull$longitude, aghy_hull$latitude, col=add.alpha('dodgerblue', 0.25),border='dodgerblue')
points(aghy$longitude,aghy$latitude, col='black', cex=0.75,pch=19)
points(centroid_aghy[,1],centroid_aghy[,2], col='red', cex=2,pch=15)

