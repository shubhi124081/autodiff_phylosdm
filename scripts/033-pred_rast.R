# rm(list = ls())
# Set up ----
HPC <- Sys.getenv("HPC")
if (HPC == FALSE) {
    root <- "~/phylo-sdms2"
} else {
    root <- "~/phylo-sdms2"
}
scripts_directory <- file.path(root, "scripts")
data_directory <- file.path(root, "data")
raw_directory <- file.path(root, "raw_data")
res_directory <- file.path(root, "res")
analysis_directory <- file.path(root, "analysis")
cond_pred_directory <- file.path(analysis_directory, "cond_pred")

# Load necessary libraries and functions
library(mvtnorm)
source(file.path(scripts_directory, "000-phyloGenie_functions.R"))
load(file.path(raw_directory, "hummingbirdSA", "spList.Rdata")) # loads obj spList
# Need SA, NA shapefile
area <- rnaturalearth::ne_countries(
    scale = "medium",
    returnclass = "sf",
    continent = c("north america", "south america")
)
# Vectorize SA and NA shapefile
area_vect <- terra::vect(area)
tempdir <- function() {
    return("/Volumes/LaCie/temp")
}
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
epath <- "~/env"

# Script will search across all clusters and subset the ones that have
# experiments defined
CLUSTERS <- names(spList)

# This will be the loop
log_file <- file.path(root, "analysis", "logs", "033-pred_rast_error_log.txt")
if (!file.exists(log_file)) {
    file.create(log_file)
}
for (i in seq_len(length(spList))) {
    cluster <- CLUSTERS[i]
    tryCatch(
        {
            print(paste0("Running cluster: ", cluster, " (", i, "/", length(spList), ")"))

            # Load raw data ----
            path <- file.path(raw_directory, "hummingbirdSA", cluster, "e2024_ALL_run_files.Rdata")
            if (file.exists(path)) {
                contents <- load(path)
            } # loads obj everything

            # Get coordinates - the LGCP runs only of presences/counts so all coordinates
            cood <- store$cood
            cood_vect <- terra::vect(cood, geom = c("lon", "lat"), crs = env_crs)

            # Get env data for the whole range
            l <- list()
            agg_fact <- 1
            widenby <- 5

            # Set range extent
            domain_ext <- terra::ext(cood_vect)
            domain_ext <- terra::extend(domain_ext, widenby)
            rm(cood_vect)

            # Env rasters
            r <- terra::rast(file.path(epath, env_files[1]))
            tmp <- terra::crop(r, domain_ext)
            tmp2 <- terra::mask(tmp, area_vect)
            r <- tmp2
            r_constant <- tmp2

            for (ii in 2:length(env_files)) {
                r1 <- terra::rast(file.path(epath, env_files[ii]))
                print(env_files[[ii]])
                tmp <- terra::crop(r1, domain_ext)
                tmp2 <- terra::mask(tmp, area_vect)
                tmp2 <- terra::resample(tmp2, r_constant)
                r <- c(r, tmp2)
            }

            # Sort out env rasters
            names(r) <- c(
                "meanTemp", "tempSeasonality", "precipQuart", "precipSeasonality",
                "cloudCover", "Annual_EVI", "TRI", "elevation"
            )

            r$meanTemp2 <- r$meanTemp^2
            r$precipQuart2 <- r$precipQuart^2

            # Create effort rasters
            x <- store$x # already scaled!
            empty <- r[[1]]
            num_observers <- empty
            terra::values(num_observers) <- median(x[, "num_observers"])
            num_observers <- terra::mask(num_observers, area_vect)
            r$num_observers <- num_observers

            duration <- empty
            terra::values(duration) <- median(x[, "duration"])
            duration <- terra::mask(duration, area_vect)
            r$duration <- duration

            distance <- empty
            terra::values(distance) <- median(x[, "distance"])
            distance <- terra::mask(distance, area_vect)
            r$distance <- distance

            Intercept <- empty
            terra::values(Intercept) <- 1
            Intercept <- terra::mask(Intercept, area_vect)
            r$Intercept <- Intercept

            # Write out the prediction raster
            terra::writeRaster(r,
                file = file.path(
                    root, "analysis",
                    "pred_rast", paste0(cluster, "_pred_rast.tiff")
                ),
                overwrite = TRUE
            )
        },
        error = function(e) {
            message <- paste0("Error in cluster: ", cluster, " - ", conditionMessage(e), "\n")
            write(message, file = log_file, append = TRUE)
            print(message)
        }
    )
}
