rm(list = ls())
# In this script, we will basically need to combine two sets of generated data as one
# write out the indices file
# write out the final data file
# There will be no reps

# Libraries and functions
source("~/phylo-sdms2/scripts/000-raw_data_gen/000-raw_data_functions.R") # nolint
source("~/phylo-sdms2/scripts/000-phyloGenie_functions.R")


# Set-up
root <- "~/phylo-sdms2"
dpath <- file.path(root, "raw_data/hummingbirdSA")
epath <- "~/env"
effpath <- file.path(dpath, "effort_rasters")

# Load clusters
load(file.path(dpath, "spList.Rdata")) # loads obj spList

# Organize
keep_info <- c("year")
keep_effort <- c("number_observers", "duration_minutes", "effort_distance_km")

# Break up by time
indices <- c(
    "TRAIN_t2003-t2018",
    "TEST_t2019-t2024"
) # Second part of the name should match the effort rasters
# Initialize a vector to store clusters with errors
error_clusters <- c()

# Loop through each cluster
#
for (cluster in names(spList)) {
    cat(sprintf("Processing cluster %s (%d out of %d)\n", cluster, which(names(spList) == cluster), length(spList)))
    tryCatch(
        {
            # Test first
            def_lev <- strsplit(indices[grep("TEST", indices)], "_")[[1]][2]
            contents <- load(file.path(dpath, cluster, paste0(def_lev, "_ALL_run_files.Rdata"))) # loads obj store
            store_test <- store # make sure to do this because store will be overwritten

            # Train second
            def_lev <- strsplit(indices[grep("TRAIN", indices)], "_")[[1]][2]
            contents <- load(file.path(dpath, cluster, paste0(def_lev, "_ALL_run_files.Rdata"))) # loads obj store

            # Number of store sites
            list_of_indices <- list()
            list_of_indices$rep_1$training_sites <- store$y$site
            list_of_indices$rep_1$testing_sites <- store_test$y$site

            save(list_of_indices, file = file.path(dpath, cluster, paste0(def_lev, "_ALL_indices.Rdata")))
        },
        error = function(e) {
            # Append the cluster name and error message to the error_clusters vector
            error_clusters <<- c(error_clusters, paste(cluster, ":", e$message))
        }
    )
}

# Print clusters with errors
if (length(error_clusters) > 0) {
    cat("The following clusters encountered errors:\n")
    cat(paste(error_clusters, collapse = "\n"), "\n")
}
