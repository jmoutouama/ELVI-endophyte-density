# remove objects and clear workspace
rm(list = ls(all=TRUE))

# Load packages
library(tidyverse)
library(terra)
library(geojsonio)
library(sp)
library(ggspatial)
library(ntbox)
library(raster)
library(rgl)
library(stringr)
library(prism)
library("RColorBrewer")

# Set seed 13
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
tmean_annual_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1993:2023))))
tmean_spring_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1993:2023, mon = 1:4))))
tmean_summer_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1993:2023, mon = 5:8))))
tmean_autumn_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1993:2023, mon = 9:12))))

# calculating standard deviation in temp
tmean_annual_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1993:2023))))
tmean_spring_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1993:2023, mon = 1:4))))
tmean_summer_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1993:2023, mon = 5:8))))
tmean_autumn_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1993:2023, mon = 9:12))))

# calculating the cumulative precipitation for each year and for each season within the year
ppt_annual <- ppt_spring <- ppt_summer <- ppt_autumn <- ppt_winter<- list()
for(y in 1993:2023){
  ppt_annual[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y))))
  ppt_spring[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y, mon = 1:4))))
  ppt_summer[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y, mon = 5:8))))
  ppt_autumn[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y, mon = 9:12)))) 
}

# Taking the mean of the cumulative precipitation values
ppt_annual_norm <- terra::mean(terra::rast(unlist(ppt_annual)))
ppt_spring_norm <- terra::mean(terra::rast(unlist(ppt_spring)))
ppt_summer_norm <- terra::mean(terra::rast(unlist(ppt_summer)))
ppt_autumn_norm <- terra::mean(terra::rast(unlist(ppt_autumn)))

# calculating the standard deviation in precip
ppt_annual_sd <- terra::stdev(terra::rast(unlist(ppt_annual)))
ppt_spring_sd <- terra::stdev(terra::rast(unlist(ppt_spring)))
ppt_summer_sd <- terra::stdev(terra::rast(unlist(ppt_summer)))
ppt_autumn_sd <- terra::stdev(terra::rast(unlist(ppt_autumn)))

# Study area shapefile 
US <- terra::vect("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/data/USA_vector_polygon/States_shapefile.shp")
US_land <- US[(!US$State_Name %in% c("HAWAII", "ALASKA", "ARIZONA", "COLORADO", "UTAH", "NEVADA", "NEW MEXICO", "IDAHO", "MONTANA", "WYOMING", "CALIFORNIA", "WASHINGTON", "OREGON")),]
US_land_reprojected <- terra::project(US_land, crs(tmean_annual_norm))
plot(US_land_reprojected)

# Crop the study area-----
tmean_spring_norm <- terra::crop(tmean_spring_norm, US_land_reprojected, mask = TRUE)
tmean_summer_norm <- terra::crop(tmean_summer_norm, US_land_reprojected, mask = TRUE)
tmean_autumn_norm <- terra::crop(tmean_autumn_norm, US_land_reprojected, mask = TRUE)
tmean_spring_sd <- terra::crop(tmean_spring_sd, US_land_reprojected, mask = TRUE)
tmean_summer_sd <- terra::crop(tmean_summer_sd, US_land_reprojected, mask = TRUE)
tmean_autumn_sd <- terra::crop(tmean_autumn_sd, US_land_reprojected, mask = TRUE)
ppt_spring_norm <- terra::crop(ppt_spring_norm, US_land_reprojected, mask = TRUE)
ppt_summer_norm <- terra::crop(ppt_summer_norm, US_land_reprojected, mask = TRUE)
ppt_autumn_norm <- terra::crop(ppt_autumn_norm, US_land_reprojected, mask = TRUE)
ppt_spring_sd <- terra::crop(ppt_spring_sd, US_land_reprojected, mask = TRUE)
ppt_summer_sd <- terra::crop(ppt_summer_sd, US_land_reprojected, mask = TRUE)
ppt_autumn_sd <- terra::crop(ppt_autumn_sd, US_land_reprojected, mask = TRUE)

## Stacking all the climatic variables 
US_land_clim <- terra::rast(list(tmean_spring_norm, tmean_summer_norm, tmean_autumn_norm,
                                 tmean_spring_sd, tmean_summer_sd, tmean_autumn_sd,
                                 ppt_spring_norm, ppt_summer_norm, ppt_autumn_norm,
                                 ppt_spring_sd, ppt_summer_sd, ppt_autumn_sd))

names(US_land_clim) <- c("tmean_spring_norm", "tmean_summer_norm", "tmean_autumn_norm",
                         "tmean_spring_sd", "tmean_summer_sd", "tmean_autumn_sd",
                         "ppt_spring_norm", "ppt_summer_norm", "ppt_autumn_norm",
                         "ppt_spring_sd", "ppt_summer_sd", "ppt_autumn_sd")

US_land_clim_stack <- stack(US_land_clim)
plot(US_land_clim_stack)

# Function to thin occurrences - Keep only one per climatic data pixel
thin_occurrences <- function(occ_data, clim_raster) {
  occ_data <- occ_data %>% 
    filter(!is.na(lat) & !is.na(lon)) %>%
    unique()
  
  # Assign raster cell ID to each occurrence point
  occ_data$cell <- terra::cellFromXY(clim_raster, occ_data[, c("lon", "lat")])
  
  # Keep only one occurrence per pixel (cell)
  occ_thinned <- occ_data %>%
    group_by(cell) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(-cell)  # Remove temporary column
  
  return(occ_thinned)
}

# Agrostis hyemalis-----
# Apply thinning function to each species dataset
aghy_occ_raw <- readRDS(url("https://www.dropbox.com/scl/fi/ijl1i7964qxzcxvm7blz8/aghy_occ_raw.rds?rlkey=msbs3xzkjc719uld8yw7hb6cj&dl=1"))
# names(aghy_occ_raw)
aghy_occ_raw %>% 
  filter(!is.na(lat) & !is.na(lon) & !is.na(year)  & as.numeric(lon >=-102.6458) & as.numeric(lat <45) & country == "United States") %>%
  unique() %>% 
  dplyr::select(country, lon, lat, year) %>% 
  arrange(lat) -> aghy

# Thin occurrences
aghy_occ_thinned <- thin_occurrences(aghy, US_land_clim)

# Plot thinned occurrence points for Agrostis hyemalis
plot(US_land_reprojected)
points(aghy_occ_thinned[, c("lon", "lat")], pch = 20, cex = 0.5, col = "red")

# Model calibration selection using Minimum Volume Ellipsoids (MVEs)
train_index_aghy <- sample(1:nrow(aghy_occ_thinned), 0.80 * nrow(aghy_occ_thinned))
test_index_aghy <- setdiff(1:nrow(aghy_occ_thinned), train_index_aghy)
# Split occurrences into train and test
aghy_train <- aghy_occ_thinned[train_index_aghy, ]
aghy_test <- aghy_occ_thinned[test_index_aghy, ]

# Extract environmental information for both train and test data
aghy_etrain <- raster::extract(US_land_clim_stack, aghy_train[, c("lon", "lat")], df = TRUE)
sum(na.omit(aghy_etrain))
aghy_etrain <- na.omit(aghy_etrain)[,-1]
summary(aghy_etrain)

aghy_etest <- raster::extract(US_land_clim_stack, aghy_test[, c("lon", "lat")], df = TRUE)
aghy_etest <- na.omit(aghy_etest)[,-1]
summary(aghy_etest)

# Find correlated environmental variables
env_varsL_aghy <- ntbox::correlation_finder(cor(aghy_etrain, method = "spearman"), threshold = 0.75, verbose = F)
env_vars_aghy <- env_varsL_aghy$descriptors
print(env_vars_aghy)

# Fit ellipsoid models
nvarstest <- c(3,4)
level <- 0.99
env_bg <- ntbox::sample_envbg(US_land_clim_stack, 20000)
omr_criteria <- 0.06
proc <- TRUE

e_select_aghy <- ntbox::ellipsoid_selection(env_train = aghy_etrain,
                                            env_test = aghy_etest,
                                            env_vars = env_vars_aghy,
                                            level = level,
                                            nvarstest = nvarstest,
                                            env_bg = env_bg,
                                            omr_criteria = omr_criteria,
                                            proc = proc)

# Display the first 10 rows of the results
head(e_select_aghy, 10)

# Best ellipsoid model for "omr_criteria"
bestvarcomb_aghy <- stringr::str_split(e_select_aghy$fitted_vars, ",")[[1]]
best_mod_aghy <- ntbox::cov_center(aghy_etrain[, bestvarcomb_aghy],
                                   mve = TRUE,
                                   level = 0.99,
                                   vars = 1:length(bestvarcomb_aghy))

# Projection model in geographic space
mProj_aghy <- ntbox::ellipsoidfit(US_land_clim_stack[[bestvarcomb_aghy]],
                                  centroid = best_mod_aghy$centroid,
                                  covar = best_mod_aghy$covariance,
                                  level = 0.95, size = 3)

if(length(bestvarcomb_aghy)==3){
  rgl::rglwidget(reuse = TRUE)
}

# Mahalanobis distance for common garden populations
garden <- read.csv("https://www.dropbox.com/scl/fi/1eu5lhkg5mx7roj3zd7g0/Study_site.csv?rlkey=tonb6sswc7zqf123ct06t64yp&dl=1", stringsAsFactors = FALSE) %>% 
  unique() %>% 
  arrange(latitude)

garden_clim <- raster::extract(US_land_clim_stack, garden[, c("longitude", "latitude")], df = TRUE)
garden_clim <- garden_clim[, -1]
mhd_aghy <- stats::mahalanobis(garden_clim[, bestvarcomb_aghy], center = best_mod_aghy$centroid, cov = best_mod_aghy$covariance)
distance_aghy <- data.frame(garden, distance = mhd_aghy)
distance_aghy$Species <- rep("AGHY", length(mhd_aghy))

plot(mhd_aghy ~ longitude, data=distance_aghy)
cor.test(distance_aghy$longitude,distance_aghy$distance)

# Elymus virginicus----
# elvi_occ_raw <- gbif(genus="Elymus",species="virginicus",download=TRUE) 
# saveRDS(elvi_occ_raw, file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/ELVI Model output/occurence/elvi_occ_raw.rds")
elvi_occ_raw <- readRDS(url("https://www.dropbox.com/scl/fi/0ssa5gepxyz28b7ykw1x8/elvi_occ_raw.rds?rlkey=4dx0q4lw2112droh73hmh7xte&dl=1"))

elvi_occ_raw %>% 
  filter(!is.na(lat) & !is.na(lon) & !is.na(year) & as.numeric(lon >= -102.6458) & country == "United States") %>% 
  unique() %>% 
  dplyr::select(country, lon, lat, year) %>% 
  arrange(lat) -> elvi

# dim(elvi)

# Thin occurrences
elvi_occ_thinned <- thin_occurrences(elvi, US_land_clim)

# Plot thinned occurrences
plot(US_land_reprojected)
points(elvi_occ_thinned[, c("lon", "lat")], pch = 20, cex = 0.5, col = "red")

# Model calibration selection using Minimum Volume Ellipsoids (MVEs).
# Random sample indexes
train_index_elvi <- sample(1:nrow(elvi_occ_thinned), 0.80 * nrow(elvi_occ_thinned))
test_index_elvi <- setdiff(1:nrow(elvi_occ_thinned), train_index_elvi)

# Split occurrences into train and test
elvi_train <- elvi_occ_thinned[train_index_elvi, ]
elvi_test <- elvi_occ_thinned[test_index_elvi, ]

# Extracts the environmental information for both train and test data
elvi_etrain <- raster::extract(US_land_clim_stack, elvi_train[, c("lon", "lat")], df = TRUE)
elvi_etrain <- na.omit(elvi_etrain)
elvi_etrain <- elvi_etrain[, -1]

elvi_etest <- raster::extract(US_land_clim_stack, elvi_test[, c("lon", "lat")], df = TRUE)
elvi_etest <- na.omit(elvi_etest)
elvi_etest <- elvi_etest[, -1]

env_varsL_elvi <- ntbox::correlation_finder(cor(elvi_etrain, method = "spearman"), threshold = 0.70, verbose = F)
env_vars_elvi <- env_varsL_elvi$descriptors
print(env_vars_elvi)

# Now we specify the number of variables to fit the ellipsoid models; in the example, we will fit for 3 dimensions
nvarstest <- c(3,4)

# Now we use the function ellipsoid_selection to run the model calibration and selection protocol
e_select_elvi <- ntbox::ellipsoid_selection(env_train = elvi_etrain,
                                            env_test = elvi_etest,
                                            env_vars = env_vars_elvi,
                                            level = level,
                                            nvarstest = nvarstest,
                                            env_bg = env_bg,
                                            omr_criteria = omr_criteria,
                                            proc = proc)

# Let’s see the first 10 rows of the results
head(e_select_elvi, 10)

# Best ellipsoid model for "omr_criteria"
bestvarcomb_elvi <- stringr::str_split(e_select_elvi$fitted_vars, ",")[[1]]

# Ellipsoid model (environmental space)
best_mod_elvi <- ntbox::cov_center(elvi_etrain[, bestvarcomb_elvi],
                                   mve = T,
                                   level = 0.99,
                                   vars = 1:length(bestvarcomb_elvi))

# Projection model in geographic space
mProj_elvi <- ntbox::ellipsoidfit(US_land_clim_stack[[bestvarcomb_elvi]],
                                  centroid = best_mod_elvi$centroid,
                                  covar = best_mod_elvi$covariance,
                                  level = 0.99, size = 3)
if(length(bestvarcomb_elvi)==3){
  rgl::rglwidget(reuse = TRUE)
}

# Mahalanobis distance for common garden populations
mhd_elvi <- stats::mahalanobis(garden_clim[, bestvarcomb_elvi], center = best_mod_elvi$centroid, cov = best_mod_elvi$covariance)
distance_elvi <- data.frame(garden, distance = mhd_elvi)
distance_elvi$Species <- rep("ELVI", length(mhd_elvi))

plot(mhd_elvi ~ longitude, data = distance_elvi)
cor.test(distance_elvi$longitude, distance_elvi$distance)

# Poa autumnalis---
# poa_occ_raw <- gbif(genus="Poa", species="autumnalis", download=TRUE) 
# saveRDS(poa_occ_raw, file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/ELVI Model output/occurence/poa_occ_raw.rds")

poa_occ_raw <- readRDS(url("https://www.dropbox.com/scl/fi/oip7ndyf0d99rqxcqxb0q/poa_occ_raw.rds?rlkey=920uql1gd4gahnh8utw9fz96l&dl=1"))

poa_occ_raw %>% 
  filter(!is.na(lat) & !is.na(lon) & !is.na(year) & as.numeric(lon >= -102.6458) & as.numeric(lat < 45) & country == "United States") %>% 
  unique() %>% 
  dplyr::select(country, lon, lat, year) %>% 
  arrange(lat) -> poa

dim(poa)

# Thin occurrences
poa_occ_thinned <- thin_occurrences(poa, US_land_clim)

# Plot thinned occurrences
plot(US_land_reprojected)
points(poa_occ_thinned[, c("lon", "lat")], pch = 20, cex = 0.5, col = "red")

# Model calibration selection using Minimum Volume Ellipsoids (MVEs).
# Random sample indexes
train_index_poa <- sample(1:nrow(poa_occ_thinned), 0.80 * nrow(poa_occ_thinned))
test_index_poa <- setdiff(1:nrow(poa_occ_thinned), train_index_poa)

# Split occurrences into train and test
poa_train <- poa_occ_thinned[train_index_poa, ]
poa_test <- poa_occ_thinned[test_index_poa, ]

# Extracts the environmental information for both train and test data
poa_etrain <- raster::extract(US_land_clim_stack, poa_train[, c("lon", "lat")], df = TRUE)
poa_etrain <- na.omit(poa_etrain)
poa_etrain <- poa_etrain[, -1]

poa_etest <- raster::extract(US_land_clim_stack, poa_test[, c("lon", "lat")], df = TRUE)
poa_etest <- na.omit(poa_etest)
poa_etest <- poa_etest[, -1]

env_varsL_poa <- ntbox::correlation_finder(cor(poa_etrain, method = "spearman"), threshold = 0.70, verbose = F)
env_vars_poa <- env_varsL_poa$descriptors
print(env_vars_poa)

# Now we specify the number of variables to fit the ellipsoid models; in the example, we will fit for 3 dimensions
nvarstest <- c(3,4)

# Now we use the function ellipsoid_selection to run the model calibration and selection protocol
e_select_poa <- ntbox::ellipsoid_selection(env_train = poa_etrain,
                                           env_test = poa_etest,
                                           env_vars = env_vars_poa,
                                           level = level,
                                           nvarstest = nvarstest,
                                           env_bg = env_bg,
                                           omr_criteria = omr_criteria,
                                           proc = proc)

# Let’s see the first 10 rows of the results
head(e_select_poa, 10)

# Best ellipsoid model for "omr_criteria"
bestvarcomb_poa <- stringr::str_split(e_select_poa$fitted_vars, ",")[[1]]

# Ellipsoid model (environmental space)
best_mod_poa <- ntbox::cov_center(poa_etrain[, bestvarcomb_poa],
                                  mve = T,
                                  level = 0.99,
                                  vars = 1:length(bestvarcomb_poa))

# Projection model in geographic space
mProj_poa <- ntbox::ellipsoidfit(US_land_clim_stack[[bestvarcomb_poa]],
                                 centroid = best_mod_poa$centroid,
                                 covar = best_mod_poa$covariance,
                                 level = 0.95, size = 3)
if(length(bestvarcomb_poa)==3){
  rgl::rglwidget(reuse = TRUE)
}
# Mahalanobis distance for common garden populations
mhd_poa <- stats::mahalanobis(garden_clim[, bestvarcomb_poa], center = best_mod_poa$centroid, cov = best_mod_poa$covariance)
distance_poa <- data.frame(garden, distance = mhd_poa)
distance_poa$Species <- rep("POAU", length(mhd_poa))

plot(mhd_poa ~ longitude, data = distance_poa)
cor.test(distance_poa$longitude, distance_poa$distance)

# Combine distance data for all species
distance_species <- bind_rows(distance_aghy, distance_elvi, distance_poa)
Species.label <- c("AGHY", "ELVI", "POAU")
names(Species.label) <- c("A. hyemalis", "E. virginicus", "P. autumnalis")


# Load necessary package
library(geosphere)

# Function to calculate the distance to the centroid using Earth's curvature
calculate_distance_to_centroid <- function(occurrence_data, garden_data) {
  
  # Filter the occurrence data (remove NA values)
  occurrence_data <- occurrence_data %>%
    filter(!is.na(lat) & !is.na(lon))
  
  # Calculate the centroid (mean latitude and longitude)
  centroid_lon <- mean(occurrence_data$lon, na.rm = TRUE)
  centroid_lat <- mean(occurrence_data$lat, na.rm = TRUE)
  
  # Calculate the distance from each garden point to the centroid using the Vincenty formula (takes curvature into account)
  garden_data$distance_to_centroid <- distVincentySphere(
    cbind(garden_data$longitude, garden_data$latitude),
    c(centroid_lon, centroid_lat)
  )
  
  return(garden_data)
}

# Example Usage:
# Assuming you have the `elvi_occ_thinned` data and `garden_clim` data

# Calculate the distance from each row in garden_clim to the centroid of elvi_occ_thinned
garden_coordinate<-garden[,1:2]
aghy_with_distance <- calculate_distance_to_centroid(aghy_occ_thinned, garden_coordinate)
plot(distance_to_centroid ~ longitude, data = aghy_with_distance)
cor.test(aghy_with_distance$longitude, aghy_with_distance$distance_to_centroid)
aghy_geo_distance<-data.frame(site_code=garden$site_code,aghy_with_distance)
aghy_geo_distance$Species <- rep("AGHY", nrow(aghy_geo_distance))
elvi_with_distance <- calculate_distance_to_centroid(elvi_occ_thinned, garden_coordinate)
plot(distance_to_centroid ~ longitude, data = elvi_with_distance)
cor.test(elvi_with_distance$longitude, elvi_with_distance$distance_to_centroid)
elvi_geo_distance<-data.frame(site_code=garden$site_code,elvi_with_distance)
elvi_geo_distance$Species <- rep("ELVI", nrow(elvi_geo_distance))
poa_with_distance <- calculate_distance_to_centroid(poa_occ_thinned, garden_coordinate)
plot(distance_to_centroid ~ longitude, data = poa_with_distance)
cor.test(poa_with_distance$longitude, poa_with_distance$distance_to_centroid)
poa_geo_distance<-data.frame(site_code=garden$site_code,poa_with_distance)
poa_geo_distance$Species <- rep("POAU", nrow(poa_geo_distance))
geo_distance<-bind_rows(aghy_geo_distance,elvi_geo_distance,poa_geo_distance)

# Load necessary packages
library(ggplot2)
library(patchwork)

# Function to create regression plots with coefficients and R2
plot_regression_with_stats <- function(data, species_name) {
  # Perform linear regression
  lm_model <- lm(distance_to_centroid ~ longitude, data = data)
  
  # Get the regression coefficient and R-squared value
  coef_val <- coef(lm_model)[2]  # Slope of the regression line
  r2_val <- summary(lm_model)$r.squared  # R-squared value
  
  # Create the plot
  p <- ggplot(data, aes(x = longitude, y = distance_to_centroid)) +
    geom_point(color = "blue", size = 1) +  # Scatter plot points
    geom_smooth(method = "lm", se = TRUE, color = "red") +  # Regression line
    labs(title = paste(species_name),
         x = "Longitude",
         y = "Distance to Centroid (m)") +
    annotate("text", x = -94, y = max(data$distance_to_centroid),
             label = paste("Slope: ", round(coef_val, 2), "\nR²: ", round(r2_val, 2)),
             hjust = 0, vjust = 1, size = 3, color = "black")
  return(p)
}

# Create plots for each species
plot_aghy <- plot_regression_with_stats(aghy_with_distance, "A. hyemalis")
plot_elvi <- plot_regression_with_stats(elvi_with_distance, "E. virginicus")
plot_poa <- plot_regression_with_stats(poa_with_distance, "P. autumnalis")

# Combine the plots into one using patchwork
combined_plot <- plot_aghy + plot_elvi + plot_poa + plot_layout(ncol = 1)

# Print the combined plot
combined_plot

distance_species<-bind_rows(distance_aghy,distance_elvi,distance_poa)
Species.label<-c("AGHY","ELVI","POAU")
names(Species.label)<-c("A. hyemalis","E. virginicus","P. autumnalis")
distance_species<-readRDS(url("https://www.dropbox.com/scl/fi/kv9j0n2pbiqgrfnm5a4wn/distance_species.rds?rlkey=vni9e8tjw9enwki0mwgnllzjc&dl=1"))

# Linear regression
# Perform correlation tests and prepare annotation data
str(distance_species)
stat_results <- distance_species %>%
  group_by(Species) %>%
  summarise(
    cor_test = list(cor.test(distance, longitude)),
    p_value = cor_test[[1]]$p.value
  ) %>%
  mutate(label = paste0("p = ", signif(p_value, 3)))

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/distance_vs_longitude_linear.pdf",useDingbats = F,height=6,width=12)
ggplot(distance_species, aes(x = longitude, y = distance))+
  labs(x="Longitude",y="Mahalanobis distance")+
  geom_point(aes(color = Species))+               
  geom_smooth(aes(color = Species, fill = Species),method="lm")+
  facet_wrap(~Species, ncol = 3, nrow = 1,labeller=labeller(Species=c("AGHY"="A. hyemalis","ELVI"="E. virginicus","POAU"="P. autumnalis")))+
  #geom_text(data = stat_results, aes(x = -94, y = 20, label = label), inherit.aes = FALSE) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  theme_bw()+
  theme(legend.position ="none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size=12, color="black",
                                    face="bold.italic"))
dev.off() 


# (Generalized Additive Models)
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/distance_vs_longitude_gam.pdf",useDingbats = F,height=6,width=12)
ggplot(distance_species, aes(x = longitude, y = distance))+
  labs(x="Longitude",y="Mahalanobis distance")+
  geom_point(aes(color = Species))+      
  geom_smooth(method = "gam", formula = y ~ s(x,k=4), aes(color = Species,fill=Species), se = TRUE,alpha = 0.2) +  # 
  facet_wrap(~Species, ncol = 3, nrow = 1,labeller=labeller(Species=c("AGHY"="A. hyemalis","ELVI"="E. virginicus","POAU"="P. autumnalis")))+
  #geom_text(data = stat_results, aes(x = -94, y = 20, label = label), inherit.aes = FALSE) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  theme_bw()+
  theme(legend.position ="none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size=12, color="black",
                                    face="bold.italic"))
dev.off()

distance_species<-cbind(distance_species,geo_distance=geo_distance[,4])
saveRDS(distance_species, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Data/distance_species.rds')


# Load required packages
library(dplyr)
library(mgcv)

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/distance_plot.pdf", width = 9, height = 8)
# Set up 2x2 layout
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2, 1))

# Custom colors and labels
cols <- c("AGHY" = "#00AFBB", "ELVI" = "#E7B800", "POAU" = "#FC4E07")
labels <- c("AGHY" = expression(italic("A. hyemalis")),
            "ELVI" = expression(italic("E. virginicus")),
            "POAU" = expression(italic("P. autumnalis")))

# Panel labels
panel_labels <- c("A", "B", "C", "D")

# First 3 plots: Longitude vs. Mahalanobis distance by species
for (i in 1:length(cols)) {
  sp <- names(cols)[i]
  dat <- subset(distance_species, Species == sp)
  
  # Fit GAM
  gam_model <- gam(distance ~ s(longitude, k = 4), data = dat)
  pred_long <- seq(min(dat$longitude), max(dat$longitude), length.out = 200)
  pred <- predict(gam_model, newdata = data.frame(longitude = pred_long), se.fit = TRUE)
  
  # Plot points
  plot(dat$longitude, dat$distance, col = cols[sp], pch = 16,
       xlab = "Longitude", ylab = "Mahalanobis distance",
       main = labels[[sp]], cex.lab = 1.3, cex.axis = 1.1)
  
  # Add smooth line with CI from GAM
  polygon(c(pred_long, rev(pred_long)),
          c(pred$fit + 2 * pred$se.fit, rev(pred$fit - 2 * pred$se.fit)),
          col = adjustcolor(cols[sp], alpha.f = 0.2), border = NA)
  lines(pred_long, pred$fit, col = cols[sp], lwd = 2)
  
  # Add panel label (A, B, C, D)
  mtext(panel_labels[i], side = 3, line = 0.2, at = par("usr")[1], adj = 0, cex = 1.25, font = 1)
}

# 4th plot: log-log plot with GAM
# Filter to valid values
valid <- distance_species$geo_distance > 0 & distance_species$distance > 0
log_geo <- log10(distance_species$geo_distance[valid])
log_dist <- log10(distance_species$distance[valid])
species_colors <- cols[distance_species$Species[valid]]

# Fit GAM to log-log
log_data <- data.frame(log_geo = log_geo, log_dist = log_dist)
gam_log <- gam(log_dist ~ s(log_geo, k = 4), data = log_data)

# Predictions
new_logx <- seq(min(log_geo), max(log_geo), length.out = 200)
pred_log <- predict(gam_log, newdata = data.frame(log_geo = new_logx), se.fit = TRUE)

# Plot log-log data
plot(log_geo, log_dist,
     col = species_colors, pch = 16,
     xlab = expression(log[10] * " Distance from geographic center"),
     ylab = expression(log[10] * " Mahalanobis distance"),
     main = "", cex.lab = 1.3, cex.axis = 1.1)

# Add confidence band and smooth from GAM
polygon(c(new_logx, rev(new_logx)),
        c(pred_log$fit + 2 * pred_log$se.fit, rev(pred_log$fit - 2 * pred_log$se.fit)),
        col = adjustcolor("black", alpha.f = 0.2), border = NA)
lines(new_logx, pred_log$fit, col = "black", lwd = 2)

# Add p-value from the GAM summary
gam_p <- summary(gam_log)$s.table[1, "p-value"]
p_label <- ifelse(gam_p < 0.001, "p < 0.001", paste0("p = ", signif(gam_p, 3)))

# Add p-value label (for GAM)
# text(x = min(log_geo) + 0.1, 
#      y = max(log_dist) - 0.1, 
#      labels = p_label, 
#      cex = 1.5, 
#      font = 1)

# Add panel label (D)
mtext(panel_labels[4], side = 3, line = 0.2, at = par("usr")[1], adj = 0, cex = 1.25, font = 1)

dev.off()


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
cuts = round(seq(0, 1, length.out=10),2)

library(viridis)

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/SDM.pdf",width=12,height=10,useDingbats = F)
par(mar=c(5,5,2,3),mfrow=c(2,2))
raster::plot(mProj_aghy$suitRaster,main="",xlab="Longitude", ylab="Latitude",cex.lab=1.5,col=viridis(99),legend=TRUE)
#points(aghy[,c("lon","lat")],pch=23,cex=0.3,col="grey")
#plot(garden_map,add=T,pch = 3,col="black",cex =2)
#plot(source_map,add=T,pch = 21,col="black",bg="red",cex =1)
mtext("A",side = 3, adj = 0,cex=1.25)
mtext(~ italic("A. hyemalis"),side = 3, adj = 0.5,cex=1.2,line=0.3)
raster::plot(mProj_elvi$suitRaster,main="",xlab="Longitude", ylab="",cex.lab=1.5,col=viridis(99),legend=TRUE)
#points(elvi[,c("lon","lat")],pch=23,cex=0.3,col="grey")
#plot(garden_map,add=T,pch = 3,col="black",cex =2)
#plot(source_map,add=T,pch = 21,col="black",bg="red",cex =1)
mtext("B",side = 3, adj = 0,cex=1.25)
mtext(~ italic("E. virginicus"),side = 3, adj = 0.5,cex=1.2,line=0.3)
raster::plot(mProj_poa$suitRaster,xlab="Longitude", ylab="Latitude",cex.lab=1.5,col=viridis(99))
#points(poa[,c("lon","lat")],pch=23,cex=0.3,col="grey")
#plot(garden_map,add=T,pch = 3,col="black",cex =2)
#plot(source_map,add=T,pch = 21,col="black",bg="red",cex =1)
mtext("C",side = 3, adj = 0,cex=1.25)
mtext(~ italic("P. autumnalis"),side = 3, adj = 0.5,cex=1.2,line=0.3)
# legend(-119, 25.5, 
#        legend=c( "GBIF occurences","Common garden sites"),
#        pch = c(23,3),
#        pt.cex=c(1.5,1.5),
#        col = c("grey50","black"),
#        pt.bg=c("grey","black"),
#        cex = 1, 
#        bty = "n", 
#        horiz = F , 
# )
dev.off()