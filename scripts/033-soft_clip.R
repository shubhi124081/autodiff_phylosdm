# Set up ----
HPC <- Sys.getenv("HPC")
if (HPC == "FALSE" || HPC == "") { # Check for both "FALSE" and empty string
    root <- "~/phylo-sdms2"
    print("Running locally")
    tempdir <- function() {
        return("/Volumes/LaCie/temp")
    }
} else {
    print("Running on HPC")
    root <- "~/phylo-sdms2"
    Sys.setenv(TMPDIR = "~/palmer_scratch")
}
print(tempdir())
DATA_SAVER <- TRUE # Big data saved on LaCie
scripts_directory <- file.path(root, "scripts")
data_directory <- file.path(root, "data")
raw_directory <- file.path(root, "raw_data")
res_directory <- file.path(root, "res")
analysis_directory <- file.path(root, "analysis")
expert_directory <- file.path(root, "expert_ranges")
if (DATA_SAVER) {
    pred_rast_directory <- file.path("/Volumes/LaCie/phylo-sdms2", "analysis", "pred_rast")
} else {
    pred_rast_directory <- file.path(analysis_directory, "pred_rast")
}

# Get all the species' name to define the exent
contents <- load(file.path(raw_directory, "hummingbirdSA", "spList.Rdata"))
source(file.path(scripts_directory, "000-phyloGenie_functions.R"))

tryCatch(
    {
        # names(spList)
        for (cluster in "Abeillia") {
            fsps <- spList[[cluster]]
            print(paste0("Soft clipping for ", cluster, " (", which(names(spList) == cluster), "/", length(spList), ")"))
            # Prediction raster
            r <- terra::rast(file.path(
                pred_rast_directory, paste0(cluster, "_", "pred_rast.tiff")
            ))

            for (i in seq_len(length(fsps))) {
                tryCatch(
                    {
                        print(paste0("Soft clipping for ", fsps[i], " (", i, "/", length(fsps), ")"))
                        fsp <- fsps[i]

                        # Get expert range ----
                        range <- terra::vect(file.path(
                            expert_directory,
                            fsp, paste0(fsp, ".shp")
                        ))
                        range_buffer <- terra::buffer(range, width = 50000) # 50km buffer
                        # 1. Create lower-res version of the raster
                        r_lowres <- terra::aggregate(r[[1]], fact = 10) # You can tweak 'fact' depending on speed/quality tradeoff

                        # 2. Calculate distance on LOW RES raster
                        distance_lowres <- terra::distance(r_lowres, range)

                        # 3. Apply logistic decay on low-res distance
                        # X0 = 100000, K = 0.0001, L = 1 - old values
                        decay_lowres <- logisticFunction(X0 = 30000, L = 1, K = 0.0002, terra::values(distance_lowres))
                        decay_lowres_raster <- terra::setValues(r_lowres, 1 - decay_lowres) # 1 - decay for soft mask

                        # 3. Resample to high-res prediction raster
                        decay_highres <- terra::resample(decay_lowres_raster, r[[1]], method = "bilinear")

                        # Rescale so that max is exactly 1
                        decay_highres <- decay_highres / terra::global(decay_highres, max, na.rm = TRUE)[1, 1]

                        # Write out rasters
                        if (DATA_SAVER) {
                            # Save to LaCie
                            output_dir <- file.path(
                                "/Volumes/LaCie/phylo-sdms2/analysis", "soft_clips"
                            )
                        } else {
                            # Save to Palmer
                            output_dir <- file.path(
                                analysis_directory, "soft_clips"
                            )
                        }
                        terra::writeRaster(
                            decay_highres,
                            file.path(output_dir, paste0(fsp, "_soft_clipped.tif")),
                            overwrite = TRUE
                        )
                    },
                    error = function(e) {
                        message(paste0("Error processing ", fsps[i], ": ", e$message))
                    }
                )
            }
        }
    },
    error = function(e) {
        message(paste0("Error processing cluster ", cluster, ": ", e$message))
    }
)
