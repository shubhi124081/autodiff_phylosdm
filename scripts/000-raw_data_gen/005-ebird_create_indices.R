# Clear workspace
rm(list = ls())

# There is a vector memory reached issue with the bigger clades
# To get around that, remove the raw dataset after it's been subsetted
# You will need to include it in the loop to read it back in for the next cluster
# Which is clunky but will work

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
TEMPORAL_DOMAIN <- "e2024" # Options: e2024, e2010, e2005, e2008
NREPS <- 5
SPLIT <- 0.75
# cluster_index <- which(names(spList) %in% c("Archilochus2"
cluster_index <- seq_along(spList) # Use all clusters
for (i in cluster_index) {
    cluster <- names(spList)[i]
    tryCatch(
        {
            contents <- load(file.path(dpath, cluster, paste0(TEMPORAL_DOMAIN, "_ALL_run_files.Rdata")))
            n <- nrow(store$x)

            # Get indices for each rep and store in a list
            list_of_indices <- list()
            for (j in seq_len(NREPS)) {
                list_of_indices[[j]] <- create_indices(store, SPLIT)
            }
            names(list_of_indices) <- paste0("rep_", 1:NREPS)
            save(list_of_indices, file = file.path(dpath, cluster, paste0(TEMPORAL_DOMAIN, "_ALL_indices.Rdata")))
        },
        error = function(e) {
            message("Error in cluster: ", cluster)
            message("Error message: ", e$message)
        }
    )
}
