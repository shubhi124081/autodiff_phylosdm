# Get the raw data generation file you want, reannotate and find mean/sd per variable and write them out

# Filepaths and directories
root <- "/Users/ss4224/phylo-sdms2"
raw_data_directory <- file.path(root, "raw_data")
dpath <- file.path(root, "raw_data/hummingbirdSA")
epath <- "~/env"

# User inputs
CLUSTER <- "Phaethornis2" # one cluster name or NULL for all
DEF_LEV <- "e2024" # Deficiency level
FSP <- "ALL" # or specific species name

# Load general functions
source(file.path(root, "scripts", "000-phyloGenie_functions.R"))
# Load species list
load(file.path(raw_data_directory, "hummingbirdSA", "spList.Rdata")) # loads spList obj

# Load cluster run files
contents <- load(file.path(raw_data_directory, "hummingbirdSA", CLUSTER, paste0(DEF_LEV, "_", FSP, "_run_files.Rdata"))) # loads store obj

# Get all coordinates
cood <- store$cood

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

rex <- annotateCoods(epath, env_files, env_crs, cood)
rex$CHELSA_bio_1_2 <- rex$CHELSA_bio_1^2
rex$CHELSA_bio_13_2 <- rex$CHELSA_bio_13^2
# Scale each column in rex and save mean and sd
rex_means <- apply(rex, 2, mean, na.rm = TRUE)
rex_sds <- apply(rex, 2, sd, na.rm = TRUE)
scales_df <- data.frame(variable = colnames(rex), mean = rex_means, sd = rex_sds)

# Repeat for effort files
reff <- annotateCoods(effpath, eff_files, env_crs, cood)
reff_means <- apply(reff, 2, mean, na.rm = TRUE)
reff_sds <- apply(reff, 2, sd, na.rm = TRUE)

# Add effort variables to scales_df
effort_df <- data.frame(
    variable = colnames(reff),
    mean = reff_means,
    sd = reff_sds
)
scales_df <- rbind(scales_df, effort_df)
scales_df$variable_name <- c(
    "meanTemp", "tempSeasonality", "precipQuart", "precipSeasonality",
    "cloudCover", "Annual_EVI", "TRI", "elevation", "meanTemp2", "precipQuart2",
    "distance", "duration", "num_observers"
)

# Save scales_df to file
write.csv(scales_df, file = file.path(raw_data_directory, "scales", paste0(CLUSTER, "_", DEF_LEV, "_", FSP, "_env_scales.csv")), row.names = FALSE)

system("say data generation complete")
