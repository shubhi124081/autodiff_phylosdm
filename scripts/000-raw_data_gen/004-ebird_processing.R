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

# Set-up
root <- "~/phylo-sdms2"
dpath <- file.path(root, "raw_data/hummingbirdSA")
epath <- "~/env"
TEMPORAL_DOMAIN <- "e2024" # Options: e2024, e2010, e2005, e2008
TEMPORAL <- FALSE # If TRUE, use temporal splits script
EXPERT_RANGE_OFFSET <- TRUE # Whether to use an expert range offset in features

# Load datasets
# Load occurrence data
load(file.path(dpath, "ebird_raw-2024.Rdata")) # loads obj raw
colnames(raw)[colnames(raw) %in% c("latitude", "longitude")] <- c("lat", "lon")
colnames(raw) <- gsub(" ", "_", colnames(raw))
raw_cood <- raw[, c("lat", "lon")]

# Complete cases only
raw <- raw[complete.cases(raw), ]

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
  "CHELSA_bio_13.tif", # precip of wettest month
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
effpath <- file.path(dpath, "effort_rasters")

# Save in a dummy variable so can recycle "tree"
tree_copy <- tree

# Organize
keep_info <- c("year")
keep_effort <- c("number_observers", "duration_minutes", "effort_distance_km")

# Break up by time
indices <- list(
  "e2024" = which(raw$year >= 2003 & raw$year <= 2024),
  "e2010" = which(raw$year >= 2003 & raw$year <= 2010),
  "e2005" = which(raw$year >= 2003 & raw$year <= 2005),
  "e2008" = which(raw$year >= 2003 & raw$year <= 2008)
)
idx <- indices[[TEMPORAL_DOMAIN]]

# Split DFs
cood_temporal <- raw[idx, c("lat", "lon")]
effort <- raw[idx, keep_effort]
info <- raw[idx, keep_info]
raw <- raw[idx, ]
save(raw, file = "~/Downloads/raw_tmp.Rdata") # random filepath

# Widen extent by
widenby <- 4
repno <- 1
# The loop with tryCatch
error_clusters <- c() # Initialize a vector to store clusters with errors
cluster_index <- which(names(spList) %in% c("Archilochus2")) # "Amazilia3", c("Archilochus2", "Phaethornis2")
# cluster_index <- seq_len(length(spList))
for (i in cluster_index) {
  cluster <- names(spList)[i]
  print(paste("Cluster", i, " : ", cluster))
  sps <- spList[[cluster]]

  tryCatch(
    {
      source("~/phylo-sdms2/scripts/000-raw_data_gen/004-ebird_processing_refactor.R")
    },
    error = function(e) {
      print(paste("Error in cluster:", cluster, " - ", e$message))
      error_clusters <- c(error_clusters, cluster) # Append cluster to error list
    }
  )
  load("~/Downloads/raw_tmp.Rdata")
}

# Print clusters with errors
if (length(error_clusters) > 0) {
  print("Clusters with errors:")
  print(error_clusters)
} else {
  print("No errors encountered.")
}

# Get rid of temp file
# system("rm ~/Downloads/raw_tmp.Rdata")

# Archive
# If you want to see stuff on a map---
# Need to load world for this
# spno <- 1
# yc_train <- cbind(cood_train, y_train)
# makeOccurrenceMap(yc_train, sps[spno])
# yc_test <- cbind(sp_test_cood[[spno]], sp_test_data[[spno]])
# makeOccurrenceMap(yc_test, sps[spno])

# tryCatch(
#   {
#     main()
#   },
#   error = function(msg) {
#     print(msg)
#     return(NA)
#   }
# )
