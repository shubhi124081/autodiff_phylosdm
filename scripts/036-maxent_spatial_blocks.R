# Start clean
rm(list = ls())

# Dependencies
library(maxnet)
library(pROC)
library(ecospat)
library(cluster)

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

# Load necessary libraries and functions
library(mvtnorm)
source(file.path(scripts_directory, "000-phyloGenie_functions.R"))
source(file.path(scripts_directory, "032-eval_functions.R"))
load(file.path(raw_directory, "hummingbirdSA", "spList.Rdata")) # loads obj spList

# Need SA, NA shapefile
area <- rnaturalearth::ne_countries(
    scale = "medium",
    returnclass = "sf",
    continent = c("north america", "south america")
)
# Vectorize SA and NA shapefile
area_vect <- terra::vect(area)
epath <- "~/env"


if (interactive()) {
    # Define the experiment
    EXP_ROOT <- "happy" # folder in analysis
    EXP_ID <- "Bkg2" # exp_root + exp_id = experiment_name
    FSP <- "ALL" # ALL if all species otherwise species name "GENUS_SPECIES"
    DEF_LEV <- "e2024" # e2024 is FULL
    REPNO <- "2" # Rep number, not number of reps
    CLUSTER <- "Abeillia" # from spList, exp_root + exp_id + cluster = shortname
    TEMPORAL <- FALSE
    TEMPORAL_TEST_DEF_LEV <- "t2009-t2024" # NULL if not temporal
} else {
    args <- commandArgs(trailingOnly = TRUE)
    EXP_ROOT <- args[1] # folder in analysis
    EXP_ID <- args[2] # exp_root + exp_id = experiment_name
    FSP <- args[3] # ALL if all species otherwise species name "GENUS_SPECIES"
    DEF_LEV <- args[4] # e2024 is FULL
    REPNO <- args[5] # Rep number, not number of reps
    CLUSTER <- args[6] # NULL if all clusters, otherwise a specific cluster
    TEMPORAL <- as.logical(args[7])
    TEMPORAL_TEST_DEF_LEV <- args[8] # NULL if not temporal
}

# Script will search across all clusters and subset the ones that have
# experiments defined
if (is.null(CLUSTER)) {
    CLUSTER <- names(spList)
}

cond_pred_directory <- file.path(analysis_directory, EXP_ROOT, "cond_pred")
eval_directory <- file.path(analysis_directory, EXP_ROOT, "eval")

# Clusters subset
results <- dir(res_directory)
focal_result_files <- results[grepl(paste0(EXP_ROOT, "_", EXP_ID), results)]
cluster_focal <- CLUSTER[sapply(CLUSTER, function(cluster) {
    any(grepl(cluster, focal_result_files))
})]

# Set tag for this run
tag <- "MAXENT"

# Load raw data
shortname <- paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER)
dataset <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO, ".Rdata")
contents <- load(file.path(data_directory, shortname, dataset)) # loads object everything

# Extract and organize
# Load the full dataset
# NOTE: CLUSTER is badly used here
contents <- load(file.path(raw_directory, "hummingbirdSA", CLUSTER, paste0(DEF_LEV, "_", FSP, "_run_files.Rdata")))

# Make sure subsetting and NA removal doesn't mess everything up
# Three matricies
# x and cood are site level, y is species-site level
y <- store$y
x <- store$x
cood <- store$cood

# Construct offset matrix
offset_vars <- c("duration", "distance", "num_observers")
x <- x[, !(colnames(x) %in% offset_vars)]

# Separate y matricies
y_species_idx <- y[, "species"]
y_values <- y[, "count"]
y_sites <- y[, "site"]
# THIN
y_values <- ifelse(y_values > 0, 1, 0)


# Species for loop starts here
# Let's say we do one species at a time
J <- length(unique(y$species))
auc_values <- numeric(J)
boyce_values <- numeric(J)
threshold <- numeric(J)
sensitivity <- numeric(J)
specificity <- numeric(J)
# === Loop over species ===
for (species_idx in seq_len(J)) {
    tryCatch(
        {
            y_species <- y[y$species == species_idx, ]
            x_species <- x[y_species$site %in% x$site, ]

            # Then we spatially block
            # Grab coordinates
            cood_train <- cood[y_species$site, ]
            coords <- cood_train[!duplicated(y_species$site), ] # Drop duplicate sites

            # Cluster into 4 spatial blocks
            set.seed(123)
            k <- 10
            spatial_blocks <- cluster::pam(coords, k)$clustering
            block_df <- data.frame(site = unique(y_species$site), block = spatial_blocks)

            # Choose one block as test set
            test_block <- c(3, 9, 1)
            test_sites <- block_df$site[block_df$block %in% test_block]
            train_sites <- block_df$site[block_df$block %in% test_block == FALSE]

            # Now we split species data into test and train
            y_train_species <- y_species[y_species$site %in% train_sites, ]
            x_train_species <- x_species[x_species$site %in% train_sites, ]
            y_test_species <- y_species[y_species$site %in% test_sites, ]
            x_test_species <- x_species[x_species$site %in% test_sites, ]

            # Split the y dataset and make binary
            y_count <- y_train_species$count
            y_bin_train <- as.numeric(ifelse(y_count > 0, 1, 0))

            # Final prep before maxent
            x_train_species <- as.data.frame(x_train_species)
            stopifnot(length(y_bin_train) == nrow(x_train_species))
            stopifnot(all(complete.cases(x_train_species)))
            stopifnot(all(!is.na(y_bin_train)))
            x_train_species <- x_train_species[, which(colnames(x_train_species) %in% c("Intercept", "meanTemp2", "precipQuart2", "site") == FALSE), drop = FALSE]
            # --- REDUCE FORMULA COMPLEXITY ---
            fml <- maxnet::maxnet.formula(
                p = y_bin_train,
                data = x_train_species,
                classes = "lq" # linear + quadratic
            )
            # --- REGULARIZED MODEL ---
            maxent_model <- maxnet::maxnet(
                p = y_bin_train,
                data = x_train_species,
                f = fml,
                regmult = 10
            )

            # Predict to test set
            x_test_species <- x_test_species[, which(colnames(x_test_species) %in% c("Intercept", "meanTemp2", "precipQuart2", "site") == FALSE), drop = FALSE]
            predicted_probs <- predict(maxent_model, newdata = x_test_species, type = "cloglog")
            predicted_probs_scaled <- predicted_probs / sum(predicted_probs, na.rm = TRUE)
            test_binary <- ifelse(y_test_species$count > 0, 1, 0)

            # AUC + threshold
            roc_obj <- pROC::roc(response = test_binary, predictor = predicted_probs_scaled)
            auc_values[species_idx] <- pROC::auc(roc_obj)

            best_thresh <- pROC::coords(roc_obj, "best", best.method = "youden", ret = "threshold")
            threshold[species_idx] <- best_thresh[[1]]

            sensitivity[species_idx] <- pROC::coords(roc_obj, x = "best", ret = "sensitivity")[[1]]
            specificity[species_idx] <- pROC::coords(roc_obj, x = "best", ret = "specificity")[[1]]

            # Boyce
            log_pred <- log(predicted_probs + 1e-10)
            log_scaled <- (log_pred - min(log_pred)) / (max(log_pred) - min(log_pred))
            obs_vals <- log(predicted_probs[test_binary == 1] + 1e-10)
            boyce_cont <- ecospat::ecospat.boyce(fit = log_scaled, obs = obs_vals, window.w = 0.1)
            boyce_values[species_idx] <- boyce_cont$cor
        },
        error = function(e) {
            cat(sprintf("Error processing species index %d: %s\n", species_idx, e$message))
            auc_values[species_idx] <- NA
            threshold[species_idx] <- NA
            boyce_values[species_idx] <- NA
            sensitivity[species_idx] <- NA
            specificity[species_idx] <- NA
        }
    )
}

# === Compile results ===
results_df <- data.frame(
    Tag = rep(tag, J),
    species_idx = 1:J,
    AUC = auc_values,
    Boyce = boyce_values,
    threshold = threshold,
    sensitivity = sensitivity,
    specificity = specificity
)
