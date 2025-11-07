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
load(file.path(raw_directory, "hummingbirdSA", "spList.Rdata")) # loads obj spList

if (interactive()) {
    # Define the experiment
    EXP_ROOT <- "happy" # folder in analysis
    EXP_ID <- "woopsie" # exp_root + exp_id = experiment_name
    FSP <- "ALL" # ALL if all species otherwise species name "GENUS_SPECIES"
    DEF_LEV <- "e2024" # e2024 is FULL
    REPNO <- "1" # Rep number, not number of reps
    CLUSTERS <- "Amazilia3" # A specific cluster or NULL
} else {
    args <- commandArgs(trailingOnly = TRUE)
    EXP_ROOT <- args[1] # folder in analysis
    EXP_ID <- args[2] # exp_root + exp_id = experiment_name
    FSP <- args[3] # ALL if all species otherwise species name "GENUS_SPECIES"
    DEF_LEV <- args[4] # e2024 is FULL
    REPNO <- args[5] # Rep number, not number of reps
    CLUSTERS <- args[6] # NULL if all clusters, otherwise a specific cluster
}
cond_pred_directory <- file.path(analysis_directory, EXP_ROOT, "cond_pred")

# Script will search across all clusters and subset the ones that have
# experiments defined
if (is.null(CLUSTERS)) {
    CLUSTERS <- names(spList)
    # Clusters subset
    results <- dir(res_directory)
    focal_result_files <- results[grepl(paste0(EXP_ROOT, "_", EXP_ID), results)]
    focal_result_files <- focal_result_files[!grepl("\\.tar\\.gz$", focal_result_files)]
    # Extract cluster names from focal_result_files
    cluster_focal <- sub(paste0("^", EXP_ROOT, "_", EXP_ID, "_"), "", focal_result_files)
} else {
    cluster_focal <- CLUSTERS
}
MODEL_NAME <- "LGCP_background" # for condpred, the model is always LGCP_background

log_directory <- file.path(analysis_directory, "logs")
if (!dir.exists(log_directory)) {
    dir.create(log_directory)
}
# I'm going to write this as a for loop but first just lay out the meat of it
# Reconstruct the name of results
for (cluster in cluster_focal) {
    error_clusters <- c() # Initialize a vector to collect clusters with errors
    tryCatch(
        {
            print(paste0("Cluster: ", cluster))
            exp_name <- paste0(EXP_ROOT, "_", EXP_ID)
            shortname <- paste0(exp_name, "_", cluster)
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

            # Some dimensions we may need
            K <- stan_data$K
            J <- stan_data$J
            N_obs <- stan_data$N_obs
            N <- stan_data$N

            # Update shortname and directories for the current cluster
            shortname <- paste0(exp_name, "_", cluster)
            data_name <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO)
            shortname_model <- paste0(data_name, "_", MODEL_NAME)
            res_dir <- file.path(res_directory, shortname)

            # Load the result file
            contents <- load(file.path(res_dir, paste0(shortname_model, ".Rdata"))) # loads obj all
            stan_fit <- all$fit
            config <- all$config

            # Load corresponding data file
            filename <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO, ".Rdata")
            contents <- load(file.path(data_directory, shortname, filename)) # loads obj everything
            stan_data <- everything$data

            # Remove site from the data
            stan_data$X <- stan_data$X[, !colnames(stan_data$X) %in% "site"]
            stan_data$K <- stan_data$K - 1

            # Some dimensions we may need
            K <- stan_data$K
            J <- stan_data$J
            N_obs <- stan_data$N_obs
            N <- stan_data$N

            # Data we need
            D_phylo <- stan_data$D_phylo # phylogenetic distance matrix
            sps <- colnames(D_phylo) # species names

            print("Starting conditional prediction...")
            B_cond_list <- list()
            for (new_index in seq_len(length(sps))) {
                observed_index <- seq_len(J)[-new_index] # Indices of species with training data
                # Call the function
                B_cond_list[[new_index]] <- cp_LGCP_simple(stan_fit, D_phylo, observed_index, new_index, stan_data)
            }
            names(B_cond_list) <- sps

            # Save the results
            filename <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO, "_cond_pred.Rdata")
            save(B_cond_list, file = file.path(cond_pred_directory, filename))
        },
        error = function(e) {
            # Print a warning and collect the cluster with error
            warning(paste("Error in cluster:", cluster, "\n", "Error message:", e$message))
            error_clusters <<- c(error_clusters, cluster)
        }
    )
}
