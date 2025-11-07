# Libraries and functions ----
source("~/phylo-sdms2/scripts/000-phyloGenie_functions.R")

# Set-up ----
dpath <- "~/phylo-sdms2/raw_data/hummingbirdSA/effort_rasters"
epath <- "~/env"

# Load datasets ----
contents <- load(file.path("~/phylo-sdms2/raw_data/hummingbirdSA/", "ebird_raw-2024.Rdata")) # loads obj raw

area <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf", continent = c("north america", "south america")
)

# Organize data -----
colnames(raw) <- gsub(" ", "_", colnames(raw))
cc <- which(complete.cases(raw) == TRUE)
raw <- raw[cc, ]

# Split effort
keep_effort <- c("number_observers", "duration_minutes", "effort_distance_km")
effort <- raw[, c("latitude", "longitude", "year", keep_effort)]
colnames(effort)[which(colnames(effort) %in% c("latitude", "longitude"))] <- c("lat", "lon")
# Split coordinates
cood <- raw[, c("latitude", "longitude")]
colnames(cood) <- c("lat", "lon")
rm(raw)

# Vectorize NA and SA shapefiles
area_vect <- terra::vect(area)

# Set up rasters to create cropped/masked env rasters
ex <- getExtentDf(cood, NAMES = c("lat", "lon"))

# Create sample raster
sample_rast <- createSampleRast(ex, area_vect, epath, "CHELSA_bio_1.tif")

# Create effort rasters
e1 <- createEffortRast(sample_rast, effort, "number_observers", PATH = dpath)
e2 <- createEffortRast(sample_rast, effort, "duration_minutes", PATH = dpath)
e3 <- createEffortRast(sample_rast, effort, "effort_distance_km", PATH = dpath)

# Create effort rasters by time
# num_observers
e1_2000 <- createEffortRast(sample_rast, effort, "number_observers",
  END_TIME = 2006, PATH = dpath
)
e1_2015 <- createEffortRast(sample_rast, effort, "number_observers",
  END_TIME = 2009, PATH = dpath
)
e1_2024 <- createEffortRast(sample_rast, effort, "number_observers",
  END_TIME = 2012, PATH = dpath
)

e1_2015 <- createEffortRast(sample_rast, effort, "number_observers",
  END_TIME = 2015, PATH = dpath
)
e1_2024 <- createEffortRast(sample_rast, effort, "number_observers",
  END_TIME = 2018, PATH = dpath
)


# duration
e2_2000 <- createEffortRast(sample_rast, effort, "duration_minutes",
  END_TIME = 2006, PATH = dpath
)
e2_2015 <- createEffortRast(sample_rast, effort, "duration_minutes",
  END_TIME = 2009, PATH = dpath
)
e2_2024 <- createEffortRast(sample_rast, effort, "duration_minutes",
  END_TIME = 2012, PATH = dpath
)

e2_2015 <- createEffortRast(sample_rast, effort, "duration_minutes",
  END_TIME = 2015, PATH = dpath
)
e2_2024 <- createEffortRast(sample_rast, effort, "duration_minutes",
  END_TIME = 2018, PATH = dpath
)

# distance
e3_2000 <- createEffortRast(sample_rast, effort, "effort_distance_km",
  END_TIME = 2006, PATH = dpath
)
e3_2015 <- createEffortRast(sample_rast, effort, "effort_distance_km",
  END_TIME = 2009, PATH = dpath
)
e3_2024 <- createEffortRast(sample_rast, effort, "effort_distance_km",
  END_TIME = 2012, PATH = dpath
)

e3_2015 <- createEffortRast(sample_rast, effort, "effort_distance_km",
  END_TIME = 2015, PATH = dpath
)
e3_2024 <- createEffortRast(sample_rast, effort, "effort_distance_km",
  END_TIME = 2018, PATH = dpath
)

# For temporal ebird splits
# distance
e3_2000 <- createEffortRast(sample_rast, effort, "effort_distance_km",
  START_TIME = 2019, END_TIME = 2024, PATH = dpath
)
# duration
e2_2000 <- createEffortRast(sample_rast, effort, "duration_minutes",
  START_TIME = 2019, END_TIME = 2024, PATH = dpath
)
# num_observers
e1_2000 <- createEffortRast(sample_rast, effort, "number_observers",
  START_TIME = 2012, END_TIME = 2024, PATH = dpath
)
# Note eX_2021 and eX are essentially the same
# These have been moved to a hummingbirdSA/effort_rasters
