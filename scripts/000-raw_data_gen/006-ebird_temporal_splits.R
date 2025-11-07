# Clear workspace
rm(list = ls())

# There is a vector memory reached issue with the bigger clades
# To get around that, remove the raw dataset after it's been subsetted
# Then read it back in for the next cluster which is clunky but works

# Libraries and functions
library(ggplot2)
library(terra)
library(rnaturalearth)
source("~/phylo-sdms2/scripts/000-raw_data_gen/000-raw_data_functions.R") # nolint
source("~/phylo-sdms2/scripts/000-phyloGenie_functions.R")
TEMPORAL <- TRUE
# Set-up
root <- "~/phylo-sdms2"
dpath <- file.path(root, "raw_data/hummingbirdSA")
epath <- "~/env"
effpath <- file.path(dpath, "effort_rasters")

# Load datasets
# Load occurrence data
load(file.path(dpath, "ebird_raw-2024.Rdata")) # loads obj raw
colnames(raw)[colnames(raw) %in% c("latitude", "longitude")] <- c("lat", "lon")
colnames(raw) <- gsub(" ", "_", colnames(raw))
raw_cood <- raw[, c("lat", "lon")]

# Complete cases only
raw <- raw[complete.cases(raw), ]
raw_copy <- raw # Make a copy of the raw dataset
save(raw_copy, file = "~/Downloads/raw_tmp.Rdata") # random filepath
rm(raw)

# Load tree
load(file.path(dpath, "mcguire_tree.Rdata")) # loads obj tree

# Load clusters
load(file.path(dpath, "spList.Rdata")) # loads obj spList

# Load shapefiles
world_vect <- terra::vect(file.path(root, "raw_data/world.shp"))
load(file.path(root, "raw_data/world.Rdata")) # loads obj world
area <- rnaturalearth::ne_countries(
    scale = "medium",
    returnclass = "sf",
    continent = c("north america", "south america")
)
# Vectorize SA and NA shapefile
area_vect <- terra::vect(area)

# Configuration
env_files <- c(
    "CHELSA_bio_1.tif", # mean annual temp
    "CHELSA_bio_4.tif", # temp seasonality
    "CHELSA_bio_13.tif", # precip of wettest 1/4
    "CHELSA_bio_15.tif", # precip seasonailty
    "cloudCover.tif", # cloud cover
    "Annual_EVI.tif", # annual evi
    "TRI.tif", # topographic ruggedness index
    "elevation_1KMmean_SRTM.tif" # elevation
)
env_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
eff_files <- c(
    "effort_rast_effort_distance_km_e2024.tiff", # distance
    "effort_rast_duration_minutes_e2024.tiff", # duration
    "effort_rast_number_observers_e2024.tiff" # num observers
)

# Save in a dummy variable so can recycle "tree"
tree_copy <- tree

# Organize
keep_info <- c("year")
keep_effort <- c("number_observers", "duration_minutes", "effort_distance_km")

# Break up by time
indices <- list(
    "TRAIN_e2003-e2018" = which(raw_copy$year >= 2003 & raw_copy$year <= 2015),
    "TEST_e2019-e2024" = which(raw_copy$year > 2003 & raw_copy$year > 2015 & raw_copy$year <= 2024)
) # Second part of the name should match the effort rasters

# TO DO: THERE IS A NAMING MISMATCH BETWEEN THIS SCRIPT AND THE TEMPORAL INDICES SCRIPT

# Widen extent by
widenby <- 4
repno <- 1
# The loop with tryCatch
error_clusters <- c() # Initialize a vector to store clusters with errors
# seq_len(length(spList))
for (i in seq_len(length(spList))) {
    cluster <- names(spList)[i]
    print(paste0(cluster))
    sps <- spList[[i]]

    for (j in seq_len(length(indices))) {
        idx <- indices[[j]]
        print(names(indices)[j])

        eff_files <- c(
            "effort_rast_effort_distance_km_e2024.tiff", # distance
            "effort_rast_duration_minutes_e2024.tiff", # duration
            "effort_rast_number_observers_e2024.tiff" # num observers
        )
        def_lev <- strsplit(names(indices)[j], "_")[[1]][2]
        eff_files <- gsub("e2024", def_lev, eff_files)
        TEMPORAL_DOMAIN <- def_lev

        tryCatch(
            {
                # Split DFs
                cood_temporal <- raw_copy[idx, c("lat", "lon")]
                effort <- raw_copy[idx, keep_effort]
                info <- raw_copy[idx, keep_info]
                raw <- raw_copy[idx, sps]
                rm(raw_copy) # Remove the raw dataset to free up memory
                source("~/phylo-sdms2/scripts/000-raw_data_gen/004-ebird_processing_refactor.R") # nolint
                load("~/Downloads/raw_tmp.Rdata")
            },
            error = function(e) {
                message(paste("Error in cluster:", cluster, "and temporal domain:", TEMPORAL_DOMAIN))
                message("Error message:", e$message)
                error_clusters <<- c(
                    error_clusters,
                    paste(cluster, TEMPORAL_DOMAIN, e$message, sep = " | ")
                )
                load("~/Downloads/raw_tmp.Rdata")
            }
        )
        load("~/Downloads/raw_tmp.Rdata")
    }
}
