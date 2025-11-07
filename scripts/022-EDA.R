# With this script we want to explore the training and
# test dataset and make sure the data going into the models
# looks good
rm(list = ls())
# Libraries ----
library(ggplot2)

# Set up ----
HPC <- Sys.getenv("HPC")
if (HPC == "FALSE") { # Check for both "FALSE" and empty string
    root <- "~/phylo-sdms2"
    print("Running locally")
} else {
    print("Running on HPC")
    root <- "~/phylo-sdms2"
}

# Directories
scripts_directory <- file.path(root, "scripts")
data_directory <- file.path(root, "data")
raw_directory <- file.path(root, "raw_data")
res_directory <- file.path(root, "res")
analysis_directory <- file.path(root, "analysis")
expert_directory <- file.path(root, "expert_ranges")

if (interactive()) {
    # Define the experiment
    EXP_ROOT <- "happy" # folder in analysis
    EXP_ID <- "new_time" # exp_root + exp_id = experiment_name
    FSP <- "ALL" # ALL if all species otherwise species name "GENUS_SPECIES"
    DEF_LEV <- "t2003-t2006" # e2024 is FULL
    REPNO <- "1" # Rep number, not number of reps
    CLUSTER <- "Aglaeactis" # A specific cluster
} else {
    args <- commandArgs(trailingOnly = TRUE)
    EXP_ROOT <- args[1] # folder in analysis
    EXP_ID <- args[2] # exp_root + exp_id = experiment_name
    FSP <- args[3] # ALL if all species otherwise species name "GENUS_SPECIES"
    DEF_LEV <- args[4] # e2024 is FULL
    REPNO <- args[5] # Rep number, not number of reps
    CLUSTER <- args[6] # A specific cluster
}
# Get all the species' name to define the exent
contents <- load(file.path(raw_directory, "hummingbirdSA", "spList.Rdata"))
source(file.path(scripts_directory, "000-phyloGenie_functions.R"))

# Area vector
area <- rnaturalearth::ne_countries(
    scale = "medium",
    returnclass = "sf",
    continent = c("north america", "south america")
)
# Vectorize SA and NA shapefile
area_vect <- terra::vect(area)

# Get training data
exp_name <- paste0(EXP_ROOT, "_", EXP_ID)
shortname <- paste0(exp_name, "_", CLUSTER)
data_name <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO)
fspALL <- spList[[CLUSTER]]

# Load all raw data
filename <- paste0(DEF_LEV, "_", FSP, "_run_files.Rdata")
contents <- load(file.path(raw_directory, "hummingbirdSA", CLUSTER, filename)) # loads obj store

# Load indicies
filename <- paste0(DEF_LEV, "_", FSP, "_indices.Rdata")
contents <- load(file.path(raw_directory, "hummingbirdSA", CLUSTER, filename)) # loads obj list_of_indices

# Look at a specific species
sps <- spList[[CLUSTER]]
plot_species_data <- function(fsps, repno) {
    # Load the expert range map for that species
    range <- terra::vect(file.path(expert_directory, fsps, paste0(fsps, ".shp")))

    # Subset the data to rep of interest
    train_sites <- list_of_indices[[repno]]$training_sites
    train_sites <- train_sites[order(train_sites)]
    test_sites <- list_of_indices[[repno]]$testing_sites
    test_sites <- test_sites[order(test_sites)]

    # Train first
    y_train <- store$y[store$y$site %in% train_sites, ]
    x <- store$x[y_train$site, , drop = FALSE]
    cood <- store$cood[y_train$site, , drop = FALSE]

    # Test second
    y_test <- store$y[store$y$site %in% test_sites, ]
    x_test <- store$x[y_test$site, , drop = FALSE]
    cood_test <- store$cood[y_test$site, , drop = FALSE]

    species_map <- store$species_index_map

    # Let's plot stuff up now
    sp_index <- as.numeric(names(species_map)[which(species_map == fsps)])
    ex <- getExtentDf(rbind(cood, cood_test), NAMES = c("lon", "lat"))
    xmin <- ex[1]
    xmax <- ex[2]
    ymin <- ex[3]
    ymax <- ex[4]

    # Set up training data
    y_train_subset <- y_train[y_train$species == sp_index, ]
    cood_train_subset <- cood[cood$site %in% y_train_subset$site, ]
    cood_train_subset <- cood_train_subset[!duplicated(cood_train_subset$site), ]
    x_train_subset <- x[x$site %in% y_train_subset$site, ]
    x_train_subset <- x_train_subset[!duplicated(x_train_subset$site), ]

    df <- cbind(cood_train_subset, y_train_subset)
    df$col <- ifelse(df[, "count"] > 0, "Presence", "Absence")
    df <- df[, c("lat", "lon", "site", "species", "count", "col")]

    # Convert range (SpatVector) to sf object
    range_sf <- sf::st_as_sf(range)
    # Plot training data
    p1 <- ggplot() +
        geom_point(
            data = df[which(df$count == 0), ], aes(x = lon, y = lat),
            fill = "#D98949", shape = 21, size = 2
        ) +
        geom_point(
            data = df[which(df$count > 0), ], aes(x = lon, y = lat),
            fill = "#7CB781", shape = 21, size = 2
        ) +
        geom_sf(
            data = area, fill = NA,
            col = "black"
        ) +
        geom_sf(
            data = range_sf, fill = NA,
            col = "red", lwd = 1.5
        ) +
        coord_sf(
            xlim = c(xmin, xmax),
            ylim = c(ymin, ymax)
        ) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        labs(fill = "") +
        ggtitle(paste0(
            "Training dataset for ", fsps, " rep no. = ", repno
        ), subtitle = paste0(
            "Number of presence points = ",
            length(df$count[df$count > 0]), "\nNumber of absence points =",
            (length(df$count[df$count == 0])
            )
        )) +
        ylab("Latitude") +
        xlab("Longitude")

    # Set up test data
    y_test_subset <- y_test[y_test$species == sp_index, ]
    cood_test_subset <- cood_test[cood_test$site %in% y_test_subset$site, ]
    cood_test_subset <- cood_test_subset[!duplicated(cood_test_subset$site), ]
    x_test_subset <- x_test[x_test$site %in% y_test_subset$site, ]
    x_test_subset <- x_test_subset[!duplicated(x_test_subset$site), ]

    df <- cbind(cood_test_subset, y_test_subset)
    df$col <- ifelse(df[, "count"] > 0, "Presence", "Absence")
    df <- df[, c("lat", "lon", "site", "species", "count", "col")]

    # Plot test data
    p2 <- ggplot() +
        geom_point(
            data = df[which(df$count == 0), ], aes(x = lon, y = lat),
            fill = "#D98949", shape = 21, size = 2
        ) +
        geom_point(
            data = df[which(df$count > 0), ], aes(x = lon, y = lat),
            fill = "#7CB781", shape = 21, size = 2
        ) +
        geom_sf(
            data = area, fill = NA,
            col = "black"
        ) +
        geom_sf(
            data = range_sf, fill = NA,
            col = "red", lwd = 1.5
        ) +
        coord_sf(
            xlim = c(xmin, xmax),
            ylim = c(ymin, ymax)
        ) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        labs(fill = "") +
        ggtitle(paste0(
            "Test dataset for ", fsps, " rep no. = ", repno
        ), subtitle = paste0(
            "Number of presence points = ",
            length(df$count[df$count > 0]), "\nNumber of absence points =",
            (length(df$count[df$count == 0])
            )
        )) +
        ylab("Latitude") +
        xlab("Longitude")

    # Return both plots
    list(training_plot = p1, test_plot = p2)
}

# Run the function for a specific species and rep number
plots <- plot_species_data(sps[1], 1)

# Display the plots
print(plots$training_plot)
print(plots$test_plot)
