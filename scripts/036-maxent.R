# Start clean
rm(list = ls())

# Dependencies
library(maxnet)
library(pROC)
library(ecospat)
library(cluster)
library(ranger)
library(PresenceAbsence)

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

# === Inputs expected ===
# store$y: long-format count data with site and species index
# store$x: covariate matrix (data.frame or matrix)
# store$cood: coordinates (not needed for MaxEnt)
# store$species_index_map: named vector mapping index -> species name
# idx_test: vector of test set site indices
# tag: string to label the evaluation set

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
    CLUSTER <- "Aglaeactis" # from spList, exp_root + exp_id + cluster = shortname
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
tag_rf <- "RF"
tag_maxent <- "MAXENT"

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
store_y <- store$y
store_x <- store$x
store_cood <- store$cood

# Load the train-test indices
contents <- load(file.path(raw_directory, "hummingbirdSA", CLUSTER, paste0(DEF_LEV, "_", FSP, "_indices.Rdata")))
idx <- list_of_indices[[as.numeric(REPNO)]]
idx_train <- idx$training_sites
idx_train <- idx_train[order(idx_train)]
# Test
idx_test <- idx$testing_sites
# idx_test <- test_sites
# idx_train <- train_sites
shared <- intersect(idx_train, idx_test)
idx_test <- setdiff(idx_test, shared)
idx_test <- idx_test[order(idx_test)]

# Subset y
y <- store_y[store_y$site %in% idx_train, ]

# Match site IDs to rows in x and cood
x <- store_x[y$site, , drop = FALSE]
cood <- store_cood[y$site, , drop = FALSE]

# Construct offset matrix
offset_vars <- c("duration", "distance", "num_observers")
offset <- x[, offset_vars]
offset <- cbind(offset, site = y$site)
x <- x[, !(colnames(x) %in% offset_vars)]
colnames(offset) <- c("duration", "distance", "num_observers", "site")

# Get ready for maxent
x_matrix <- x

# Separate y matricies
y_species_idx <- y[, "species"]
y_values <- y[, "count"]
y_sites <- y[, "site"]
# THIN
y_values <- ifelse(y_values > 0, 1, 0)

# Others
offset_matrix <- offset
offset_multiplied <- offset_matrix[, "duration"] * offset_matrix[, "distance"] * offset_matrix[, "num_observers"]
offset_multiplied <- cbind(offset_multiplied, offset_matrix[, "site"])
colnames(offset_multiplied) <- c("offset", "site")
offset_vector <- as.vector(offset_multiplied[which(offset_multiplied[, "site"] %in% y_sites), "offset"])

# I don't know what's going on here - compress
y <- store$y
x <- store$x
cood <- store$cood
species_index_map <- store$species_index_map
# Remove site
if ("site" %in% colnames(x)) {
    x <- as.matrix(x[, colnames(x) != "site"]) # Shape: (N_test, K)
}
x <- x[, !(colnames(x) %in% offset_vars)]

# Separate test data
y_test <- y[y$site %in% idx_test, ]
x_test <- x[y_test$site, , drop = FALSE]
cood_test <- cood[y_test$site, , drop = FALSE]

# Process offset and remove offset vars
# offset_vars <- c("duration", "distance", "num_observers")
# processed_offset <- process_offset(x_test, y_test, offset_vars)
# offset_test <- processed_offset$offset_test
# x_test <- processed_offset$x_test
if ("site" %in% colnames(x_test)) {
    x_test <- as.matrix(x_test[, colnames(x_test) != "site"]) # Shape: (N_test, K)
}

# Prepare results containers
J <- length(unique(y_test$species))
auc_values_maxent <- numeric(J)
boyce_values_maxent <- numeric(J)
threshold_maxent <- numeric(J)
sensitivity_maxent <- numeric(J)
specificity_maxent <- numeric(J)

auc_values_rf <- numeric(J)
boyce_values_rf <- numeric(J)
threshold_rf <- numeric(J)
sensitivity_rf <- numeric(J)
specificity_rf <- numeric(J)

# === Loop over species ===
for (species_idx in seq_len(J)) {
    tryCatch(
        {
            cat("Processing species index:", species_idx, "\n")

            # Filter train + test data
            y_test_species <- y_test[y_test$species == species_idx, ]
            site_idx_test <- y_test_species$site
            x_test_species <- x[site_idx_test, , drop = FALSE]
            x_test_species <- as.data.frame(x_test_species[, !(colnames(x_test_species) %in% offset_vars)])
            x_test_species <- x_test_species[, c(
                "meanTemp", "tempSeasonality",
                "precipQuart", "precipSeasonality", "elevation",
                "Annual_EVI", "TRI", "cloudCover"
            )]

            y_train_species <- y[y$species == species_idx, ]
            site_idx_train <- y_train_species$site
            x_train_species <- x[site_idx_train, , drop = FALSE]
            y_bin_train <- ifelse(y_train_species$count > 0, 1, 0)
            x_train_species <- as.data.frame(x_train_species[, !(colnames(x_train_species) %in% offset_vars)])
            x_train_species <- x_train_species[, c(
                "meanTemp", "tempSeasonality",
                "precipQuart", "precipSeasonality", "elevation",
                "Annual_EVI", "TRI", "cloudCover"
            )]

            # Sanity check for minimum presence
            if (sum(y_bin_train) < 1) {
                warning("Fewer than 1 presences — skipping species.")
                next
            }
            y_bin_train <- as.numeric(y_train_species$count > 0)
            # Final prep before maxent
            x_train_species <- as.data.frame(x_train_species)
            stopifnot(length(y_bin_train) == nrow(x_train_species))
            stopifnot(all(complete.cases(x_train_species)))
            stopifnot(all(!is.na(y_bin_train)))
            x_train_species <- x_train_species[, colnames(x_train_species) != "Intercept", drop = FALSE]
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
                regmult = 1
            )

            pred_names <- colnames(x_train_species)
            train_data <- cbind(x_train_species, y_bin_train)
            colnames(train_data) <- c(pred_names, "pres_abs")
            # formula
            m_formula <- paste("pres_abs ~", paste(pred_names, collapse = "+"))
            # model
            rf_occ <- ranger::ranger(
                formula = m_formula,
                num.trees = 100,
                importance = "impurity",
                # num.threads = n_threads,
                respect.unordered.factors = "order",
                always.split.variables = NULL,
                probability = TRUE,
                replace = TRUE,
                data = train_data
            )

            # Predict to test set
            # random forest
            x_test_species <- x_test_species[, colnames(x_test_species) != "Intercept", drop = FALSE]
            predicted_probs_rf <- predict(rf_occ, data = x_test_species, type = "response")$predictions[, 2]
            predicted_probs_scaled_rf <- predicted_probs_rf / sum(predicted_probs_rf, na.rm = TRUE)

            # maxent
            predicted_probs_maxent <- predict(maxent_model, newdata = x_test_species, type = "cloglog")
            predicted_probs_scaled <- predicted_probs_maxent / sum(predicted_probs_maxent, na.rm = TRUE)
            test_binary <- ifelse(y_test_species$count > 0, 1, 0)

            # AUC + threshold
            roc_obj <- pROC::roc(response = test_binary, predictor = predicted_probs_scaled)
            auc_values_maxent[species_idx] <- pROC::auc(roc_obj)
            roc_obj_rf <- pROC::roc(response = test_binary, predictor = predicted_probs_scaled_rf)
            auc_values_rf[species_idx] <- pROC::auc(roc_obj_rf)

            best_thresh <- pROC::coords(roc_obj, "best", best.method = "youden", ret = "threshold")
            threshold_maxent[species_idx] <- best_thresh[[1]]

            best_thresh <- pROC::coords(roc_obj_rf, "best", best.method = "youden", ret = "threshold")
            threshold_rf[species_idx] <- best_thresh[[1]]

            sensitivity_maxent[species_idx] <- pROC::coords(roc_obj, x = "best", ret = "sensitivity")[[1]]
            specificity_maxent[species_idx] <- pROC::coords(roc_obj, x = "best", ret = "specificity")[[1]]

            sensitivity_rf[species_idx] <- pROC::coords(roc_obj_rf, x = "best", ret = "sensitivity")[[1]]
            specificity_rf[species_idx] <- pROC::coords(roc_obj_rf, x = "best", ret = "specificity")[[1]]

            # Boyce
            log_pred <- log(predicted_probs_maxent + 1e-10)
            log_scaled <- (log_pred - min(log_pred)) / (max(log_pred) - min(log_pred))
            obs_vals <- log(predicted_probs_maxent[test_binary == 1] + 1e-10)
            boyce_cont <- ecospat::ecospat.boyce(fit = log_scaled, obs = obs_vals, window.w = 0.1)
            boyce_values_maxent[species_idx] <- boyce_cont$cor

            log_pred <- log(predicted_probs_rf + 1e-10)
            log_scaled <- (log_pred - min(log_pred)) / (max(log_pred) - min(log_pred))
            obs_vals <- log(predicted_probs_rf[test_binary == 1] + 1e-10)
            boyce_cont <- ecospat::ecospat.boyce(fit = log_scaled, obs = obs_vals, window.w = 0.1)
            boyce_values_rf[species_idx] <- boyce_cont$cor
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
results_df_maxent <- data.frame(
    Tag = rep(tag_maxent, J),
    species_idx = 1:J,
    AUC = auc_values_maxent,
    Boyce = boyce_values_maxent,
    threshold = threshold_maxent,
    sensitivity = sensitivity_maxent,
    specificity = specificity_maxent
)

results_df_rf <- data.frame(
    Tag = rep(tag_rf, J),
    species_idx = 1:J,
    AUC = auc_values_rf,
    Boyce = boyce_values_rf,
    threshold = threshold_rf,
    sensitivity = sensitivity_rf,
    specificity = specificity_rf
)

# Map species names
results_df$species <- species_index_map[as.character(results_df$species_idx)]

# === Save output ===
eval_filename <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO, "_maxent_eval.Rdata")
save(results_df, file = file.path(eval_directory, eval_filename))

# Optional: print and inspect
print(results_df)

r <- terra::rast("~/phylo-sdms2/analysis/pred_rast/Aglaeactis_pred_rast.tiff")
r <- r[[c(
    colnames(x_train_species)
)]]

predict_fun <- function(df) {
    df <- as.data.frame(df)
    df <- df[, colnames(x_train_species), drop = FALSE] # Ensure order
    predict(maxent_model, newdata = df, type = "cloglog")
}

all_na_layers <- names(r)[sapply(1:terra::nlyr(r), function(i) all(is.na(terra::values(r[[i]]))))]
if (length(all_na_layers) > 0) warning("Some layers have only NA: ", paste(all_na_layers, collapse = ", "))

prediction_raster <- terra::predict(
    object = r,
    model = maxent_model,
    type = "cloglog",
    na.rm = TRUE
)
