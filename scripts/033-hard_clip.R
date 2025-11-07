# Spatial pred directory
spatial_directory <- file.path("/vast/palmer/scratch/jetz/ss4224", "spatial_pred")
# spatial_directory <- file.path("/Volumes/LaCie/phylo-sdms2", "analysis", "spatial_pred")
expert_directory <- file.path("/vast/palmer/pi/jetz/ss4224/phylo-sdms2/expert_ranges")
# expert_directory <- file.path("~/phylo-sdms2/expert_ranges")

# Specs of the files we need
EXP_ROOT <- "happy"
EXP_ID <- "runfast15"
DEF_LEV <- "t2015"
REPNO <- "1"
EXTRA <- "_local_"
CORRECTED <- "CORRECTED" # or NOT_CORRECTED

# List all files
all_files <- list.files(spatial_directory, full.names = TRUE)

# Filter based on all patterns
# Build a regex pattern to match the required filename structure
pattern <- paste0(
    "^", EXP_ROOT, "_", EXP_ID, "_.*_", DEF_LEV, "_", REPNO, "_.*_", CORRECTED, EXTRA, "\\.tif$"
)
matched_files <- all_files[grepl(pattern, basename(all_files))]

# Loop through matched files
# Get name of the species
# File structure is EXP_ROOT_EXP_ID_CLUSTER_FSP_DEF_LEV_REPNO_EXTRA.tif
# Get species name from the file name
for (file in matched_files) {
    # Extract the species name from the file name
    parts <- strsplit(basename(file), "_")[[1]]
    cluster <- parts[3]
    fsp <- parts[c(4, 5)]
    fsp <- paste(fsp, collapse = "_") # Join the parts to get the full species name
    tag <- paste(parts[8:(length(parts) - 1)], collapse = "_")

    # Try-catch block for processing each file
    tryCatch(
        {
            # Load the raster
            r <- terra::rast(file)

            # Get the associated range
            range_file <- file.path(expert_directory, fsp, paste0(fsp, ".shp"))
            if (file.exists(range_file)) {
                range <- terra::vect(range_file)

                # Hard clip the raster using the range
                r_clipped <- terra::mask(r, range)

                # Save the clipped raster
                output_file <- file.path(spatial_directory, paste0(
                    EXP_ROOT, "_", EXP_ID, "_", cluster, "_",
                    fsp, "_", DEF_LEV, "_", REPNO, "_", tag, "_hard_clipped.tif"
                ))
                terra::writeRaster(r_clipped, output_file, overwrite = TRUE)

                print(paste("Saved hard clipped raster for", fsp, "to", output_file))
            } else {
                print(paste("Range file not found for", fsp))
            }
        },
        error = function(e) {
            print(paste("Error processing", fsp, ":", e$message))
        }
    )
}
