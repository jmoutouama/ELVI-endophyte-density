#### Project: 
#### PURPOSE: Estimating distance from climatic niche center (Mahalanobis distance)
#### AUTHOR: Jacob Moutouama
#### DATE LAST MODIFIED: 
# remove objects and clear workspace
rm(list = ls(all=TRUE))
#load packages
# devtools::install_github("luismurao/ntbox")
library(tidyverse)
library(terra)
library(geojsonio)
library(sp)
library(rgeos)
library(ggspatial)
library(ntbox)
library(raster)
library(rgl)
library(stringr)
library(prism)
library("RColorBrewer")
#Set seed 13
set.seed(13)
# Climatic data----
## Data from PRISM
# making a folder to store prism data
options(prism.path = '/Users/jm200/Documents/Prism Range limit/')
# getting monthly data for mean temp and precipitation
# takes a long time the first time, but can skip when you have raster files saved on your computer.
# get_prism_monthlys(type = "tmean", years = 1990:2024, mon = 1:12, keepZip = FALSE)
# get_prism_monthlys(type = "ppt", years = 1990:2024, mon = 1:12, keepZip = FALSE)
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

# Study area shapefile 
US<-terra::vect("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/data/USA_vector_polygon/States_shapefile.shp")
plot(US)
US_land<-US[(!US$State_Name %in% c("HAWAII","ALASKA")),]
plot(US_land)

# Crop the study area-----
tmean_spring_norm<-terra::crop(tmean_spring_norm, US_land,mask=TRUE)
tmean_summer_norm<-terra::crop(tmean_summer_norm, US_land,mask=TRUE)
tmean_autumn_norm<-terra::crop(tmean_autumn_norm, US_land,mask=TRUE)
tmean_spring_sd<-terra::crop(tmean_spring_sd, US_land,mask=TRUE)
tmean_summer_sd<-terra::crop(tmean_summer_sd, US_land,mask=TRUE)
tmean_autumn_sd<-terra::crop(tmean_autumn_sd, US_land,mask=TRUE)
ppt_spring_norm<-terra::crop(ppt_spring_norm, US_land,mask=TRUE)
ppt_summer_norm<-terra::crop(ppt_summer_norm, US_land,mask=TRUE)
ppt_autumn_norm<-terra::crop(ppt_autumn_norm, US_land,mask=TRUE)
ppt_spring_sd<-terra::crop(ppt_spring_sd, US_land,mask=TRUE)
ppt_summer_sd<-terra::crop(ppt_summer_sd, US_land,mask=TRUE)
ppt_autumn_sd<-terra::crop(ppt_autumn_sd, US_land,mask=TRUE)

##  Stacking all the climatic variables 
US_land_clim<-terra::rast(list(tmean_spring_norm,tmean_summer_norm,tmean_autumn_norm,
                                    tmean_spring_sd,tmean_summer_sd,tmean_autumn_sd,
                                    ppt_spring_norm,ppt_summer_norm,ppt_autumn_norm,
                                    ppt_spring_sd,ppt_summer_sd,ppt_autumn_sd))
names(US_land_clim) <- c("tmean_spring_norm","tmean_summer_norm","tmean_autumn_norm",
                              "tmean_spring_sd","tmean_summer_sd","tmean_autumn_sd",
                              "ppt_spring_norm","ppt_summer_norm","ppt_autumn_norm",
                              "ppt_spring_sd","ppt_summer_sd","ppt_autumn_sd")
US_land_clim_stack<-stack(US_land_clim)
plot(US_land_clim_stack)

# Agrostis hyemalis-----
#aghy_occ_raw <- gbif(genus="Agrostis",species="hyemalis",download=TRUE) 
aghy_occ_raw<-readRDS(url("https://www.dropbox.com/scl/fi/ijl1i7964qxzcxvm7blz8/aghy_occ_raw.rds?rlkey=msbs3xzkjc719uld8yw7hb6cj&dl=1"))
names(aghy_occ_raw) 
aghy_occ_raw %>% 
  filter(!is.na(lat) & !is.na(lon) & !is.na(year) & year %in% (1950:2024) & country=="United States" ) %>%
  unique() %>% 
  dplyr::select(country,lon, lat,year)%>% 
  arrange(lat)->aghy
#dim(aghy)
plot(US_land)
points(aghy[,c("lon","lat")],pch=20,cex=0.5,col="red")
# Model calibration selection using Minimum Volume Ellipsoids (MVEs).
# Random sample indexes
train_index_aghy <- sample(1:nrow(aghy), 0.80 * nrow(aghy))
test_index_aghy <- setdiff(1:nrow(aghy), train_index_aghy)
## Split occurences in train and test
aghy_train <- aghy[train_index_aghy, ]
aghy_test <- aghy[test_index_aghy, ]

# Extracts the environmental information for both train and test data
aghy_etrain <- raster::extract(US_land_clim_stack,aghy_train[,c("lon", "lat")],df=TRUE)
aghy_etrain<-na.omit(aghy_etrain)
aghy_etrain <- aghy_etrain[,-1]
head(aghy_etrain)

aghy_etest <- raster::extract(US_land_clim_stack,aghy_test[,c("lon","lat")], df=TRUE)
aghy_etest<-na.omit(aghy_etest)
aghy_etest <- aghy_etest[,-1]
head(aghy_etest)

env_varsL_aghy <- ntbox::correlation_finder(cor(aghy_etrain, method = "spearman"),threshold = 0.70,verbose = F)
env_vars_aghy <- env_varsL_aghy$descriptors
print(env_vars_aghy )

#Now we specify the number of variables to fit the ellipsoid models; in the example, we will fit for 3,5, and 6 dimensions
nvarstest <- 3
## This parameter is to specify the proportion of training points that will be used to fit the minimum volume ellipsoid (Van Aelst and Rousseeuw 2009).
# Level
level <- 0.99
# This background data is just to compute the partial ROC test
env_bg <- ntbox::sample_envbg(US_land_clim_stack,10000)
## For selecting the model we will use an arbitrary value of 6 percent of omission; it is not a rule but accepted omission rates are those bellow 10%. We will ask the function to return the partial ROC value (Peterson, Papes, and Soberon 2008)
omr_criteria <- 0.05
proc <- TRUE

# Now we just need to use the function ellipsoid_selection to run the model calibration and selection protocol
e_select_aghy <- ntbox::ellipsoid_selection(env_train = aghy_etrain,
                                            env_test = aghy_etest,
                                            env_vars = env_vars_aghy,
                                            level = level,
                                            nvarstest = nvarstest,
                                            env_bg = env_bg,
                                            omr_criteria= omr_criteria,
                                            proc = proc)


# Let’s see the first 10 rows of the results
head(e_select_aghy,10)
# With the following lines of code, I am going to display the model in the first row of the table
# Best ellipsoid model for "omr_criteria" 
bestvarcomb_aghy <- stringr::str_split(e_select_aghy$fitted_vars,",")[[1]]
# Ellipsoid model (environmental space)
best_mod_aghy <- ntbox::cov_center(aghy_etrain[,bestvarcomb_aghy],
                                   mve = T,
                                   level = 0.99,
                                   vars = 1:length(bestvarcomb_aghy))


# Projection model in geographic space
#install.packages("rgl",dependencies = TRUE)
mProj_aghy <- ntbox::ellipsoidfit(US_land_clim_stack[[bestvarcomb_aghy]],
                             centroid = best_mod_aghy$centroid,
                             covar = best_mod_aghy$covariance,
                             level = 0.99,size = 3)

# Elymus virginicus---
#elvi_occ_raw <- gbif(genus="Elymus",species="virginicus",download=TRUE) 
#saveRDS(elvi_occ_raw, file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/ELVI Model output/occurence/elvi_occ_raw.rds")
elvi_occ_raw<-readRDS(url("https://www.dropbox.com/scl/fi/0ssa5gepxyz28b7ykw1x8/elvi_occ_raw.rds?rlkey=4dx0q4lw2112droh73hmh7xte&dl=1"))
elvi_occ_raw %>% 
  filter(!is.na(lat) & !is.na(lon) & !is.na(year) & year %in% (1950:2024) & country=="United States") %>% 
  unique() %>% 
  dplyr::select(country,lon, lat,year)%>% 
  arrange(lat)->elvi
dim(elvi)
plot(US_land)
points(elvi[,c("lon","lat")],pch=20,cex=0.1,col="red")
# Model calibration selection using Minimum Volume Ellipsoids (MVEs).
# Random sample indexes
train_index_elvi <- sample(1:nrow(elvi), 0.80 * nrow(elvi))
test_index_elvi <- setdiff(1:nrow(elvi), train_index_elvi)
# # Split occurences in train and test
elvi_train <- elvi[train_index_elvi, ]
elvi_test <- elvi[test_index_elvi, ]

# Extracts the environmental information for both train and test data
elvi_etrain <- raster::extract(US_land_clim_stack,elvi_train[,c("lon", "lat")],df=TRUE)
elvi_etrain<-na.omit(elvi_etrain)
elvi_etrain <- elvi_etrain[,-1]
head(elvi_etrain)

elvi_etest <- raster::extract(US_land_clim_stack,elvi_test[,c("lon","lat")], df=TRUE)
elvi_etest<-na.omit(elvi_etest)
elvi_etest <- elvi_etest[,-1]
head(elvi_etest)

env_varsL_elvi <- ntbox::correlation_finder(cor(elvi_etrain, method = "spearman"),threshold = 0.70,verbose = F)
env_vars_elvi <- env_varsL_elvi$descriptors
print(env_vars_elvi )

#Now we specify the number of variables to fit the ellipsoid models; in the example, we will fit for 3,5, and 6 dimensions
nvarstest <- c(3,4)

# Now we just need to use the function ellipsoid_selection to run the model calibration and selection protocol
e_select_elvi <- ntbox::ellipsoid_selection(env_train = elvi_etrain,
                                            env_test = elvi_etest,
                                            env_vars = env_vars_elvi,
                                            level = level,
                                            nvarstest = nvarstest,
                                            env_bg = env_bg,
                                            omr_criteria= omr_criteria,
                                            proc = proc)


# Let’s see the first 10 rows of the results
head(e_select_elvi,10)
# With the following lines of code, I am going to display the model in the first row of the table
# Best ellipsoid model for "omr_criteria" 
bestvarcomb_elvi <- stringr::str_split(e_select_elvi$fitted_vars,",")[[1]]
# Ellipsoid model (environmental space)
best_mod_elvi <- ntbox::cov_center(elvi_etrain[,bestvarcomb_elvi],
                                   mve = T,
                                   level = 0.99,
                                   vars = 1:length(bestvarcomb_elvi))


# Projection model in geographic space
mProj_elvi <- ntbox::ellipsoidfit(US_land_clim_stack[[bestvarcomb_elvi]],
                             centroid = best_mod_elvi$centroid,
                             covar = best_mod_elvi$covariance,
                             level = 0.99,size = 3)

# Poa autumnalis---
#poa_occ_raw <- gbif(genus="Poa",species="autumnalis",download=TRUE) 
#saveRDS(poa_occ_raw, file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/ELVI Model output/occurence/poa_occ_raw.rds")
poa_occ_raw<-readRDS(url("https://www.dropbox.com/scl/fi/oip7ndyf0d99rqxcqxb0q/poa_occ_raw.rds?rlkey=920uql1gd4gahnh8utw9fz96l&dl=1"))
poa_occ_raw %>% 
  filter(!is.na(lat) & !is.na(lon) & !is.na(year) & year %in% (1950:2024) & country=="United States") %>% 
  unique() %>% 
  dplyr::select(country,lon, lat,year)%>% 
  arrange(lat)->poa
dim(poa)
# plot(US_land)
# points(poa[,c("lon","lat")],pch=20,cex=0.1)

# Model calibration selection using Minimum Volume Ellipsoids (MVEs).
# Random sample indexes
train_index_poa <- sample(1:nrow(poa), 0.70 * nrow(poa))
test_index_poa <- setdiff(1:nrow(poa), train_index_poa)
# # Split occurences in train and test
poa_train <- poa[train_index_poa, ]
poa_test <- poa[test_index_poa, ]

# Extracts the environmental information for both train and test data
poa_etrain <- raster::extract(US_land_clim_stack,poa_train[,c("lon", "lat")],df=TRUE)
poa_etrain<-na.omit(poa_etrain)
poa_etrain <- poa_etrain[,-1]
head(poa_etrain)

poa_etest <- raster::extract(US_land_clim_stack,poa_test[,c("lon","lat")], df=TRUE)
poa_etest<-na.omit(poa_etest)
poa_etest <- poa_etest[,-1]
head(poa_etest)

env_varsL_poa <- ntbox::correlation_finder(cor(poa_etrain, method = "spearman"),threshold = 0.70,verbose = F)
env_vars_poa <- env_varsL_poa$descriptors
print(env_vars_poa )

#Now we specify the number of variables to fit the ellipsoid models; in the example, we will fit for 3,5, and 6 dimensions
nvarstest <- 3

# Now we just need to use the function ellipsoid_selection to run the model calibration and selection protocol
e_select_poa <- ntbox::ellipsoid_selection(env_train = poa_etrain,
                                            env_test = poa_etest,
                                            env_vars = env_vars_poa,
                                            level = level,
                                            nvarstest = nvarstest,
                                            env_bg = env_bg,
                                            omr_criteria= omr_criteria,
                                            proc = proc)


# Let’s see the first 10 rows of the results
head(e_select_poa,10)
# With the following lines of code, I am going to display the model in the first row of the table
# Best ellipsoid model for "omr_criteria" 
bestvarcomb_poa <- stringr::str_split(e_select_poa$fitted_vars,",")[[1]]
# Ellipsoid model (environmental space)
best_mod_poa <- ntbox::cov_center(poa_etrain[,bestvarcomb_poa],
                                   mve = T,
                                   level = 0.99,
                                   vars = 1:length(bestvarcomb_poa))


# Projection model in geographic space")
mProj_poa <- ntbox::ellipsoidfit(US_land_clim_stack[[bestvarcomb_poa]],
                             centroid = best_mod_poa$centroid,
                             covar = best_mod_poa$covariance,
                             level = 0.99,size = 3)

# mahalanobis distance 
read.csv("https://www.dropbox.com/scl/fi/1eu5lhkg5mx7roj3zd7g0/Study_site.csv?rlkey=tonb6sswc7zqf123ct06t64yp&dl=1", stringsAsFactors = F) %>% 
  unique() %>% 
  arrange(latitude)->garden ## common garden populations

garden_clim<- raster::extract(US_land_clim_stack,garden[,c("longitude","latitude")], df=TRUE)
garden_clim <- garden_clim[,-1]
mhd_aghy <- stats::mahalanobis(garden_clim[,bestvarcomb_aghy],center = best_mod_aghy$centroid,cov = best_mod_aghy$covariance)
distance_aghy<-data.frame(garden,distance=mhd_aghy)
distance_aghy$Species<-rep("AGHY",length(mhd_aghy))

plot(mhd_aghy ~ longitude, data=distance_aghy)
cor.test(distance_aghy$longitude,distance_aghy$distance)

mhd_elvi <- stats::mahalanobis(garden_clim[,bestvarcomb_elvi],center = best_mod_elvi$centroid,cov = best_mod_elvi$covariance)
distance_elvi<-data.frame(garden,distance=mhd_elvi)
distance_elvi$Species<-rep("ELVI",length(mhd_elvi))

plot(mhd_elvi ~ longitude, data=distance_elvi)
cor.test(distance_elvi$longitude,distance_elvi$distance)

mhd_poa <- stats::mahalanobis(garden_clim[,bestvarcomb_poa],center = best_mod_poa$centroid,cov = best_mod_poa$covariance)
distance_poa<-data.frame(garden,distance=mhd_poa)
distance_poa$Species<-rep("POAU",length(mhd_poa))

plot(mhd_poa ~ longitude, data=distance_poa)
cor.test(distance_poa$longitude,distance_poa$distance)

distance_species<-bind_rows(distance_aghy,distance_elvi,distance_poa)
Species.label<-c("AGHY","ELVI","POAU")
names(Species.label)<-c("A. hyemalis","E. virginicus","P. autumnalis")

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/distance_vs_longitude.pdf",useDingbats = F,height=6,width=12)
ggplot(distance_species, aes(x = longitude, y = distance))+
  labs(x="Longitude",y="Mahalanobis distance")+
  geom_point(aes(color = Species))+               
  geom_smooth(aes(color = Species, fill = Species))+
  facet_wrap(~Species, ncol = 3, nrow = 1,labeller=labeller(Species=c("AGHY"="A. hyemalis","ELVI"="E. virginicus","POAU"="P. autumnalis")))+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  theme_bw()+
  theme(legend.position ="none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
    strip.text.x = element_text(size=12, color="black",
                                   face="bold.italic"))
dev.off() 

saveRDS(distance_species, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Data/distance_species.rds')


read.csv("https://www.dropbox.com/scl/fi/1eu5lhkg5mx7roj3zd7g0/Study_site.csv?rlkey=tonb6sswc7zqf123ct06t64yp&dl=1", stringsAsFactors = F) %>% 
  unique() %>% 
  arrange(latitude)->garden_map ## common garden population
read.csv("https://www.dropbox.com/scl/fi/go448bqe9z6meisgkd6g9/source_pop.csv?rlkey=oeutzzf4lgo616zeam06tul8t&dl=1", stringsAsFactors = F) %>% 
  dplyr::select(latitude,longitude) %>% 
  unique() %>% 
  arrange(latitude)->source_map ## source populations

sp::coordinates(garden_map) <- ~ longitude + latitude
sp::coordinates(source_map) <- ~ longitude + latitude
CRS1 <- CRS("+init=epsg:4326") # WGS 84
crs(garden_map) <- CRS1
crs(source_map) <- CRS1
cuts = round(seq(0, 1, length.out=20),2)


pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/SDM.pdf",width=12,height=10,useDingbats = F)
par(mar=c(5,5,2,3),mfrow=c(2,2))
raster::plot(mProj_aghy$suitRaster,main="",xlab="Longitude", ylab="Latitude",cex.lab=1.5,breaks=cuts,col=terrain.colors(20),legend=FALSE)
points(aghy[,c("lon","lat")],pch=23,cex=0.3,col="grey")
plot(garden_map,add=T,pch = 3,col="black",cex =2)
#plot(source_map,add=T,pch = 21,col="black",bg="red",cex =1)
mtext("A",side = 3, adj = 0,cex=1.25)
mtext(~ italic("A. hyemalis"),side = 3, adj = 0.5,cex=1.2,line=0.3)
raster::plot(mProj_elvi$suitRaster,main="",xlab="Longitude", ylab="",cex.lab=1.5,breaks=cuts,col=terrain.colors(20),legend=FALSE)
points(elvi[,c("lon","lat")],pch=23,cex=0.3,col="grey")
plot(garden_map,add=T,pch = 3,col="black",cex =2)
#plot(source_map,add=T,pch = 21,col="black",bg="red",cex =1)
mtext("B",side = 3, adj = 0,cex=1.25)
mtext(~ italic("E. virginicus"),side = 3, adj = 0.5,cex=1.2,line=0.3)
raster::plot(mProj_poa$suitRaster,xlab="Longitude", ylab="Latitude",cex.lab=1.5,breaks=cuts,col=terrain.colors(20))
points(poa[,c("lon","lat")],pch=23,cex=0.3,col="grey")
plot(garden_map,add=T,pch = 3,col="black",cex =2)
#plot(source_map,add=T,pch = 21,col="black",bg="red",cex =1)
mtext("C",side = 3, adj = 0,cex=1.25)
mtext(~ italic("P. autumnalis"),side = 3, adj = 0.5,cex=1.2,line=0.3)
legend(-119, 25.5, 
       legend=c( "GBIF occurences","Common garden sites"),
       pch = c(23,3),
       pt.cex=c(1.5,1.5),
       col = c("grey50","black"),
       pt.bg=c("grey","black"),
       cex = 1, 
       bty = "n", 
       horiz = F , 
)
dev.off()