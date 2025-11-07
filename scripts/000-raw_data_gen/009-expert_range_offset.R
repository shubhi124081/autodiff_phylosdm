APPROACH <- 1 #  1 to be sourced by 004-ebird_processing_refactor.R when EXPERT_RANGE_OFFSET is TRUE
# 2 # Instead of regenerating a whole dataset, we can load an existing one and just modify x

if (APPROACH == 1) {
    # We need to go species by species and exp(soft_clip) values and store them as an additional column in x
    soft_clip_list <- list()
    for (fsps in sps) {
        tryCatch(
            {
                raster_path <- file.path(offset_path, paste0(fsps, "_soft_clipped.tif"))

                # Get species index and sites first (used whether raster exists or not)
                sp_index <- which(species_index_map == fsps)
                y_count_species <- y_count[y_count$species == sp_index, , drop = FALSE]
                y_species_sites <- y_count_species$site
                cood_species <- cood[cood$site %in% y_species_sites, , drop = FALSE]

                if (!file.exists(raster_path)) {
                    # If raster missing, fill soft_clip with 1 for all species' sites
                    warning(sprintf("soft clip raster not found for '%s'; filling with 1s: %s", fsps, raster_path))
                    if (nrow(cood_species) > 0) {
                        soft_clip_list[[fsps]] <- data.frame(site = cood_species$site, soft_clip = rep(1, nrow(cood_species)))
                    } else {
                        soft_clip_list[[fsps]] <- data.frame(site = integer(0), soft_clip = numeric(1))
                    }
                    next
                }

                # Load in species soft clip
                soft_clip_rast <- terra::rast(raster_path)

                if (nrow(cood_species) == 0) stop("no coordinates for species sites")

                # Annotate soft clip values
                ext <- terra::extract(soft_clip_rast, cood_species[, c("lon", "lat")], method = "bilinear")
                if (is.null(ext) || ncol(ext) < 2) stop("unexpected extract result")
                soft_clip_values <- ext[, 2]

                # Build a data frame
                soft_clip_df <- data.frame(
                    site = cood_species$site,
                    soft_clip = soft_clip_values
                )
                soft_clip_list[[fsps]] <- soft_clip_df
            },
            error = function(e) {
                warning(sprintf("Failed processing species '%s': %s", fsps, e$message))
                # fallback: create NA rows for sites (if any) so merge won't drop them
                sp_index <- which(species_index_map == fsps)
                y_count_species <- y_count[y_count$species == sp_index, , drop = FALSE]
                y_species_sites <- y_count_species$site
                cood_species <- cood[cood$site %in% y_species_sites, , drop = FALSE]
                if (nrow(cood_species) > 0) {
                    soft_clip_list[[fsps]] <- data.frame(site = cood_species$site, soft_clip = NA_real_)
                } else {
                    soft_clip_list[[fsps]] <- data.frame(site = integer(0), soft_clip = numeric(0))
                }
                invisible(NULL)
            }
        )
    }

    # Now need to merge soft_clip values into x
    tryCatch(
        {
            soft_clip_all <- do.call(rbind, soft_clip_list)
            if (nrow(soft_clip_all) > 0) {
                soft_clip_all <- aggregate(soft_clip ~ site, data = soft_clip_all, FUN = max)
                soft_clip_all$soft_clip <- soft_clip_all$soft_clip + 0.00001 # small constant to avoid zeros
                x <- merge(x, soft_clip_all, by = "site", all.x = TRUE)
            } else {
                warning("No soft_clip values available; adding soft_clip column with NA")
                x$soft_clip <- 1
            }
        },
        error = function(e) {
            warning("Failed to merge soft_clip values into x: ", e$message)
            x$soft_clip <- 1
        }
    )
}

if (APPROACH == 2) {
    # Libraries and functions
    library(ggplot2)
    library(terra)
    library(rnaturalearth)
    source("~/phylo-sdms2/scripts/000-raw_data_gen/000-raw_data_functions.R") # nolint
    source("~/phylo-sdms2/scripts/000-phyloGenie_functions.R")

    offset_path <- "/Volumes/LaCie/phylo-sdms2/analysis/soft_clips"
    if (!dir.exists(offset_path)) {
        warning("Expert range offset path not found. Please ensure Volumes/LaCie is connected.")
    }

    # Set-up
    root <- "~/phylo-sdms2"
    dpath <- file.path(root, "raw_data/hummingbirdSA")
    epath <- "~/env"
    TEMPORAL_DOMAIN <- "e2024" # Options: e2024, e2010, e2005, e2008
    FSP <- "ALL"

    # Load clusters
    load(file.path(dpath, "spList.Rdata")) # loads obj spList
    CLUSTERS <- "Colibri" # one cluster or names(spList)

    for (CLUSTER in CLUSTERS) {
        raw_file <- file.path(
            dpath, CLUSTER, paste0(DEF_LEV, "_", FSP, "_run_files.Rdata")
        )
        contents <- load(raw_file)
        species_map <- store$species_index_map
        y <- store$y
        x <- store$x
        cood <- store$cood

        # Now we can add the soft_clip values to x
        sps <- store$tree$tip.label
        soft_clip_all <- data.frame()
        soft_clip_list <- list()
        for (fsps in sps) {
            cat("Processing species:", fsps, "\n")
            # Load in species soft clip
            soft_clip_rast <- terra::rast(file.path(offset_path, paste0(fsps, "_soft_clipped.tif")))

            # Get species' sites
            fsps_index <- as.numeric(names(species_map)[which(species_map == fsps)])
            y_species <- y[y$species == fsps_index, ]
            y_species_sites <- y_species$site

            # Get coordinates for the sites
            cood_species <- cood[cood$site %in% y_species_sites, , drop = FALSE]

            # Annotate soft clip values
            soft_clip_values <- terra::extract(soft_clip_rast, cood_species[, c("lon", "lat")], method = "bilinear")[, 2]

            # Build a data frame
            soft_clip_df <- data.frame(
                site = cood_species$site,
                soft_clip = soft_clip_values + 0.00001 # small constant to avoid zeros
            )
            soft_clip_list[[fsps]] <- soft_clip_df
        }

        soft_clip_all <- do.call(rbind, soft_clip_list)

        # Merge into x using site
        x <- merge(x, soft_clip_all, by = "site", all.x = TRUE)

        # Save modified store
        store$x <- x
        save(store, file = raw_file)
    }
}
