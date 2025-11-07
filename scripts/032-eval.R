rm(list = ls())

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

if (interactive()) {
    # Define the experiment
    EXP_ROOT <- "happy" # folder in analysis
    EXP_ID <- "woopsie" # exp_root + exp_id = experiment_name
    FSP <- "ALL" # ALL if all species otherwise species name "GENUS_SPECIES"
    DEF_LEV <- "e2024" # e2024 is FULL
    REPNO <- "1" # Rep number, not number of reps
    CLUSTERS <- "Archilochus2" # from spList, exp_root + exp_id + cluster = shortname
    TEMPORAL <- FALSE # TRUE if temporal, FALSE if not
    TEMPORAL_TEST_DEF_LEV <- "t2019-t2024" # NULL if not temporal
    ART <- FALSE # Artificial thinning experiment
} else {
    args <- commandArgs(trailingOnly = TRUE)
    EXP_ROOT <- args[1] # folder in analysis
    EXP_ID <- args[2] # exp_root + exp_id = experiment_name
    FSP <- args[3] # ALL if all species otherwise species name "GENUS_SPECIES"
    DEF_LEV <- args[4] # e2024 is FULL
    REPNO <- args[5] # Rep number, not number of reps
    CLUSTERS <- args[6] # NULL if all clusters, otherwise a specific cluster
    TEMPORAL <- as.logical(args[7])
    TEMPORAL_TEST_DEF_LEV <- args[8] # NULL if not temporal
    ART <- as.logical(args[9]) # Artificial thinning experiment
}

set.seed(10212025)

# Script will search across all clusters and subset the ones that have
# experiments defined
if (is.null(CLUSTERS)) {
    CLUSTERS <- names(spList)
}

cond_pred_directory <- file.path(analysis_directory, EXP_ROOT, "cond_pred")
eval_directory <- file.path(analysis_directory, EXP_ROOT, "eval")

# Clusters subset
results <- dir(res_directory)
focal_result_files <- results[grepl(paste0(EXP_ROOT, "_", EXP_ID), results)]
cluster_focal <- CLUSTERS[sapply(CLUSTERS, function(cluster) {
    any(grepl(cluster, focal_result_files))
})]

# TAGS LOOP STARTS HERE
TAGS <- c("PHYLO_BASE", "PHYLO_ALL", "NOPHYLO_BASE")
for (tag in TAGS) {
    if (tag == "NOPHYLO_BASE") {
        MODEL_NAME <- "LGCP_nophylo_background"
    } else {
        MODEL_NAME <- "LGCP_background"
    }

    # Reconstruct the name of results
    exp_name <- paste0(EXP_ROOT, "_", EXP_ID)
    shortname <- paste0(exp_name, "_", cluster_focal)
    data_name <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO)

    shortname_model <- paste0(data_name, "_", MODEL_NAME)
    res_dir <- file.path(res_directory, shortname)
    # Load the result file
    contents <- load(file.path(res_dir, paste0(shortname_model, ".Rdata"))) # loads obj all
    stan_fit <- all$fit
    config <- all$config

    # Load corresponding data file
    contents <- load(file.path(data_directory, shortname, paste0(data_name, ".Rdata"))) # loads obj everything
    stan_data <- everything$data

    # Load the full dataset
    if (TEMPORAL) {
        contents <- load(file.path(
            raw_directory, "hummingbirdSA",
            cluster_focal,
            paste0(TEMPORAL_TEST_DEF_LEV, "_", FSP, "_run_files.Rdata")
        )) # loads obj store
    } else {
        contents <- load(file.path(
            raw_directory, "hummingbirdSA", cluster_focal,
            paste0(DEF_LEV, "_", FSP, "_run_files.Rdata")
        )) # loads obj store
    }

    # Need the test data indices
    contents <- load(file.path(raw_directory, "hummingbirdSA", cluster_focal, paste0(DEF_LEV, "_", FSP, "_indices.Rdata"))) # loads list_of_indices
    idx <- list_of_indices[[as.numeric(REPNO)]]
    idx_test <- idx$testing_sites
    idx_test <- idx_test[order(idx_test)]

    # Create test data ---------
    # Y
    y <- store$y
    y_test <- y[y$site %in% idx_test, ]

    # X
    x <- store$x
    x_test <- x[y_test$site, ]
    if ("site" %in% colnames(x_test)) {
        x_test <- as.matrix(x_test[, colnames(x_test) != "site"]) # Shape: (N_test, K)
    }

    # Cood
    cood <- store$cood
    cood_test <- cood[y_test$site, ]

    # Offset
    # if ("soft_clip" %in% colnames(x_test)) {
    #     # offset_vars <- c("duration", "distance", "num_observers", "soft_clip")
    #     offset_vars <- "soft_clip"
    # } else {
    #     offset_vars <- c("duration", "distance", "num_observers")
    # }
    # offset_vars <- "soft_clip
    offset_vars <- config$data$offsets
    processed_offset <- process_offset(x_test, y_test, offset_vars)
    offset_test <- processed_offset$offset_test
    offset_test <- offset_test * 0
    x_test <- processed_offset$x_test

    # Dimensions needed
    N_test <- nrow(x_test) # Number of test locations
    K <- ncol(x_test) # Number of environmental predictors
    J <- length(unique(y_test$species))

    # Extract posterior samples
    posterior_B <- rstan::extract(stan_fit, pars = "B")$B # (num_samples, J, K)
    posterior_sigma_f <- rstan::extract(stan_fit, pars = "sigma_f")$sigma_f # (num_samples)

    if (tag == "PHYLO_ALL") {
        cond_filename <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO, "_cond_pred.Rdata")
        contents <- load(file.path(cond_pred_directory, cond_filename)) # loads obj B_cond_list
        # contents <- load("/Users/ss4224/phylo-sdms2/analysis/cond_pred/happy_test_PBkg_Aglaeactis_ALL_e2024_1_cond_pred.Rdata")
        # Convert B_cond_list to an array of dimensions (550, 3, 11)
        num_samples <- dim(posterior_B)[1] # 550
        num_species <- dim(posterior_B)[2] # 3
        num_predictors <- dim(posterior_B)[3] # 11

        # Create an empty array to hold the combined values
        B_cond_array <- array(NA, dim = c(num_samples, num_species, num_predictors))

        # Fill the array using a loop over species
        if (ART) {
            species_indices <- seq_len(length(B_cond_list))
        } else {
            species_indices <- seq_len(J)
        }
        for (j in species_indices) {
            B_cond_array[, j, ] <- as.matrix(B_cond_list[[j]])
        }

        posterior_B <- B_cond_array # (num_samples, J, K)
    }
    if (ART) {
        species_indices <- unique(y_test$species) # Use species indices from y_test
    } else {
        species_indices <- 1:J
    }

    # This is where the for loop for the species index is needed
    # Predicted intensities at  locations
    auc_values <- numeric(length(species_indices))
    boyce_values <- numeric(length(species_indices))
    threshold <- numeric(length(species_indices))
    sensitivity <- numeric(length(species_indices))
    specificity <- numeric(length(species_indices))
    for (ii in seq_len(length(species_indices))) {
        species_idx <- species_indices[ii]
        tryCatch(
            {
                y_test_species <- y_test[which(y_test$species == species_idx), ]
                x_test_species <- x_test[which(y_test$species == species_idx), ]
                offset_test_species <- offset_test[which(y_test$species == species_idx)]
                # x_test_species <- x_test_species[, c("Intercept", "meanTemp", "meanTemp2", "tempSeasonality", "precipQuart",
                # "precipQuart2", "precipSeasonality", "elevation", "Annual_EVI", "TRI", "cloudCover")]

                predicted_intensity <- predict_intensity_rf(posterior_B, posterior_sigma_f, x_test_species, y_test_species, offset_test_species)
                predicted_intensity_scaled <- predicted_intensity / sum(predicted_intensity, na.rm = TRUE)

                test_binary <- ifelse(y_test$count[y_test$species == species_idx] > 0, 1, 0)

                # Calculate AUC + associated metrics
                roc_obj <- pROC::roc(response = test_binary, predictor = predicted_intensity_scaled)
                auc_values[ii] <- pROC::auc(roc_obj)
                best_thresh <- pROC::coords(roc_obj, "best", best.method = "youden", ret = "threshold")
                threshold[ii] <- sapply(best_thresh, function(x) x[1])
                sens <- pROC::coords(roc_obj, x = "best", ret = "sensitivity")
                sensitivity[ii] <- sapply(sens, function(x) x[1])
                specs <- pROC::coords(roc_obj, x = "best", ret = "specificity")
                specificity[ii] <- sapply(specs, function(x) x[1])

                # Calculate Boyce Index
                log_pred <- log(predicted_intensity + 1e-10) # Avoid -Inf
                log_scaled <- (log_pred - min(log_pred)) / (max(log_pred) - min(log_pred))
                # boyce_cont <- ecospat::ecospat.boyce(fit = log_scaled, obs = log(predicted_intensity[y_test_species$count > 0] + 1e-10), window.w = 0.1)
                # boyce_values[ii] <- boyce_cont$cor
                boyce_values[ii] <- NA
            },
            error = function(e) {
                message(sprintf("Error processing species index %d: %s", species_idx, e$message))
                auc_values[ii] <- NA
                threshold[ii] <- NA
                boyce_values[ii] <- NA
            }
        )
    }

    # Append results to a data frame
    new_results <- data.frame(
        Tag = rep(tag, length(species_indices)),
        species_idx = species_indices,
        AUC = auc_values,
        Boyce = boyce_values,
        threshold = threshold,
        sensitivity = as.numeric(sensitivity),
        specificity = as.numeric(specificity)
    )

    if (exists("results_df")) {
        results_df <- rbind(results_df, new_results)
    } else {
        results_df <- new_results
    }
}

# Get species index map
contents <- load(file.path(
    raw_directory, "hummingbirdSA", cluster_focal,
    paste0(DEF_LEV, "_", FSP, "_run_files.Rdata")
)) # loads obj store
species_idx_map <- store$species_index_map
results_df$species <- species_idx_map[match(results_df$species_idx, names(species_idx_map))]

# Save the results
eval_filename <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO, "_eval", ".Rdata")
save(results_df, file = file.path(eval_directory, eval_filename))
print(results_df)
# Archive
# EDA
# library(ggplot2)
# load(file.path(root, "raw_data/world.Rdata")) # loads obj world
# species_idx <- 2
# y_test_species <- y_test[which(y_test$species == species_idx), ]
# y_test_species_binary <- y_test_species
# y_test_species_binary$count <- ifelse(y_test_species_binary$count > 0, 1, 0)
# cood_test_species <- cood_test[which(y_test$species == species_idx), ]
# yc_test_species <- cbind(cood_test_species, y_test_species_binary)
# colnames(yc_test_species) <- c("y", "x", "site", "species", "count")
# makeOccurrenceMap(yc_test_species, FSP = "count")

# Other archived stuff
# cood_pres <- cood_test
# Y_plot <- rep(c(1, 0), c(nrow(cood_pres), nrow(cood_bg)))
# YC <- cbind(rbind(cood_pres, cood_bg), Y_plot)
# colnames(YC) <- c("lat", "lon", "A")
# makeOccurrenceMap(YC, "A", area)

# # Functions
# predict_intensity_from_stan <- function(stan_fit, x_test, y_test, offset) {
#     # Extract posterior samples
#     posterior_B <- rstan::extract(stan_fit, pars = "B")$B # (num_samples, J, K)
#     posterior_f <- rstan::extract(stan_fit, pars = "f")$f # (num_samples, N_obs)

#     num_samples <- dim(posterior_B)[1]
#     N_test <- nrow(x_test)

#     log_lambda_samples <- matrix(NA, nrow = num_samples, ncol = N_test)
#     pb_samples <- txtProgressBar(min = 0, max = num_samples, style = 3)

#     for (i in 1:num_samples) {
#         for (n in 1:N_test) {
#             species_n <- y_test$species[n] # species index
#             log_lambda_samples[i, n] <- sum(x_test[n, ] * posterior_B[i, species_n, ]) + posterior_f[i, n] + offset[n]
#         }
#         setTxtProgressBar(pb_samples, i)
#     }

#     close(pb_samples)

#     # Compute mean predicted intensity (lambda)
#     # predicted_intensity <- exp(rowMeans(log_lambda_samples))

#     return(predicted_intensity)
# }

# process_x_true_test_2 <- function(x_true_test, cood_test = NULL) {
#     x_true_test[is.infinite(x_true_test[, "distance"]), "distance"] <- 0
#     # Step 1: Add column names
#     colnames(x_true_test) <- c(
#         "meanTemp", "tempSeasonality", "precipQuart",
#         "precipSeasonality", "cloudCover", "Annual_EVI", "TRI", "elevation",
#         "distance", "duration", "num_observers"
#     )

#     # Step 2: Add new predictors
#     x_true_test$meanTemp2 <- (x_true_test$meanTemp)^2
#     x_true_test$precipQuart2 <- (x_true_test$precipQuart)^2

#     # Step 3: Scale the data
#     x_true_test <- scale_safe(x_true_test)

#     # Step 4: Convert to data frame and add Intercept
#     x_true_test <- as.data.frame(x_true_test)
#     x_true_test$Intercept <- rep(1, nrow(x_true_test))

#     # Step 5: Convert back to matrix
#     x_true_test <- as.matrix(x_true_test)

#     # Step 6: Reorder columns
#     order <- c(
#         "Intercept", "meanTemp", "meanTemp2", "tempSeasonality",
#         "precipQuart", "precipQuart2", "precipSeasonality", "elevation",
#         "Annual_EVI", "TRI", "cloudCover", "duration", "distance",
#         "num_observers"
#     )
#     x_true_test <- x_true_test[, order]

#     # Step 7: Remove rows with missing values
#     rm_ind <- which(complete.cases(x_true_test) == FALSE)
#     if (length(rm_ind) > 1) {
#         x_true_test <- x_true_test[-rm_ind, ]
#         if (!is.null(cood_test)) {
#             cood_test <- cood_test[-rm_ind, ]
#         }
#     }

#     # Return the processed data
#     if (is.null(cood_test)) {
#         return(x_true_test)
#     } else {
#         return(list(x_true_test = x_true_test, cood_test = cood_test))
#     }
# }
# thresholds <- seq(0, 0.1, length.out = 100)
# sens_spec <- sapply(thresholds, function(thresh) {
#     pred_binary <- ifelse(predicted_intensity > thresh, 1, 0)
#     cm <- table(factor(test_binary, levels = c(0, 1)), factor(pred_binary, levels = c(0, 1)))
#     sens <- cm[2, 2] / (cm[2, 1] + cm[2, 2]) # Sensitivity
#     spec <- cm[1, 1] / (cm[1, 1] + cm[1, 2]) # Specificity
#     return(c(sens = sens, spec = spec))
#     # spec <- cm[1, 1] / (cm[1, 1] + cm[1, 2])
#     # return(c(sens = sens, spec = spec))
# })
# # Find the threshold that maximizes the sum of sensitivity and specificity
# best_thresh <- thresholds[which.max(colSums(sens_spec))]
