# Set up ----
HPC <- Sys.getenv("HPC")
root <- ifelse(HPC == FALSE, "~/phylo-sdms2", "~/project/phylo-sdms2")
raw_dir <- file.path(root, "raw_data/hummingbirdSA")
expert_dir <- file.path(root, "expert_ranges")

# Libraries and functions ----
# 000-phyloGenie.R will hold all the phyloGenie functions
source(file.path(root, "scripts/000-phyloGenie_functions.R"))
require(ggplot2)

# User-input specifications ----
SP <- "Phaethornis mexicanus" # GENUS [SPACE] SPECIES
CLUSTER <- "Chlorostilbon" # Cluster it belongs to
DATASET <- "ebird_raw-2024.Rdata"
ENV_CRS <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Load in data ----
contents <- load(file.path(raw_dir, DATASET)) # loads obj raw

# Process begins here ---
# Grab species occ column with coods
cood_species <- raw[, c("latitude", "longitude", SP)]
colnames(cood_species) <- c("latitude", "longitude", "sp")
cood_species1 <- cood_species[which(cood_species$sp == 1), ]

# Make vector ----
cood_species_vect <- terra::vect(cood_species1, geom = c("longitude", "latitude"))
convex_hull <- terra::convHull(cood_species_vect)
terra::crs(convex_hull) <- ENV_CRS
# Increase the range by 50 km ----
convex_hull2 <- terra::buffer(convex_hull, width = 50000) # 50 km buffer

# Plot to make sure if you want ----
terra::plot(convex_hull, col = "blue")
terra::points(cood_species_vect)
# terra::lines(area_vect, col = "red")
terra::lines(convex_hull2, col = "blue")


# Save ----
SP <- gsub(" ", "_", SP)
terra::writeVector(convex_hull2, paste0("expert_ranges", SP, "/", paste0(SP, ".shp")))
