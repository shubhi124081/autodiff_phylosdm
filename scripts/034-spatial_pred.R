rm(list = ls())

# ---- Setup ----
HPC <- Sys.getenv("HPC")
DATA_SAVER <- TRUE
EXP_ROOT <- "happy"

if (HPC == "FALSE" || HPC == "") {
    root <- "~/phylo-sdms2"
    print("Running locally")
    tempdir <- function() "/Volumes/LaCie/temp"
    output_dir <- file.path("/Volumes/LaCie/phylo-sdms2", "analysis", EXP_ROOT, "spatial_pred")
} else {
    print("Running on HPC")
    root <- "/vast/palmer/pi/jetz/ss4224/phylo-sdms2"
    output_dir <- file.path("/vast/palmer/scratch/jetz/ss4224/spatial_pred")
    tempdir <- function() "/vast/palmer/scratch/jetz/ss4224"
}

if (dir.exists("/Volumes/LaCie/temp")) {
    terra::terraOptions(tempdir = "/Volumes/LaCie/temp")
} else {
    message("External temp drive not found. Using default tempdir.")
}

print(tempdir())
print(output_dir)

# ---- Utility Functions ----
cleanup_terra_temp <- function(path = terra::terraOptions()$tempdir, days_old = 0) {
    files <- list.files(path, full.names = TRUE)
    old_files <- files[file.info(files)$mtime < Sys.time() - days_old * 86400]
    unlink(old_files, recursive = TRUE)
    message(length(old_files), " files removed from terra tempdir.")
}

# ---- Directory Paths ----
scripts_directory <- file.path(root, "scripts")
data_directory <- file.path(root, "data")
raw_directory <- file.path(root, "raw_data")
res_directory <- file.path(root, "res")
analysis_directory <- file.path(root, "analysis")
expert_directory <- file.path(root, "expert_ranges")

if (DATA_SAVER) {
    pred_rast_directory <- if (HPC == "FALSE" || HPC == "") {
        file.path("/Volumes/LaCie/phylo-sdms2", "analysis", "pred_rast")
    } else {
        file.path("/vast/palmer/pi/jetz/ss4224/phylo-sdms2/analysis/pred_rast")
    }
} else {
    pred_rast_directory <- file.path(analysis_directory, "pred_rast")
}

# ---- Experiment Parameters ----
# Set experiment parameters either interactively or via command line arguments

if (interactive()) {
    # If running interactively (e.g., in RStudio), set parameters manually
    EXP_ROOT <- "happy" # Root name for experiment
    EXP_ID <- "clean_alloffset" # Experiment ID
    FSP <- "ALL" # Focal species or "ALL"
    DEF_LEV <- "e2024" # Definition level
    REPNO <- "1" # Replicate number
    CLUSTER <- "Colibri" # Cluster name
    HARD_CLIP <- FALSE # Whether to hard clip predictions to expert range
    BINARY_CLIP <- FALSE # Whether to output binary hard-clipped raster
    BINARY_NONCLIP <- FALSE # Whether to output binary non-clipped raster
    SOFT_CLIP <- FALSE # Whether to output soft-clipped raster
    SOFT_CLIP_IN_PRED <- TRUE # Whether to soft-clip predictions in model training
    WIDEN_BY <- 2 # How much to widen the raster plot beyond the range box
} else {
    # If running via command line, parse arguments
    args <- commandArgs(trailingOnly = TRUE)
    EXP_ROOT <- args[1] # Root name for experiment
    EXP_ID <- args[2] # Experiment ID
    FSP <- args[3] # Focal species or "ALL"
    DEF_LEV <- args[4] # Definition level
    REPNO <- args[5] # Replicate number
    CLUSTER <- args[6] # Cluster name
    HARD_CLIP <- as.logical(args[7]) # Whether to hard clip predictions to expert range
    BINARY_CLIP <- as.logical(args[8]) # Whether to output binary hard-clipped raster
    BINARY_NONCLIP <- as.logical(args[9]) # Whether to output binary non-clipped raster
    SOFT_CLIP <- as.logical(args[10])
    SOFT_CLIP_IN_PRED <- as.logical(args[11]) # Whether to soft-clip predictions in model training
    WIDEN_BY <- as.numeric(args[12]) # How much to widen the raster
}

# ---- Load Data ----
load(file.path(raw_directory, "hummingbirdSA", "spList.Rdata"))
source(file.path(scripts_directory, "000-phyloGenie_functions.R"))

exp_name <- paste0(EXP_ROOT, "_", EXP_ID)
shortname <- paste0(exp_name, "_", CLUSTER)
data_name <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO)
fspALL <- spList[[CLUSTER]]
cond_pred_directory <- file.path(analysis_directory, EXP_ROOT, "cond_pred")

# ---- Prediction Raster ----
r <- terra::rast(file.path(pred_rast_directory, paste0(CLUSTER, "_pred_rast.tiff")))

# ---- Environmental Variables ----
env_vars <- c(
    "Intercept", "meanTemp", "meanTemp2", "tempSeasonality",
    "precipQuart", "precipQuart2", "precipSeasonality", "elevation",
    "Annual_EVI", "TRI", "cloudCover"
)
r_env <- r[[env_vars]]
K <- dim(r_env)[3]

# Standardize environmental layers
# Load in means/sds
scales_df <- read.csv(file.path(raw_directory, "hummingbirdSA", CLUSTER, paste0(DEF_LEV, "_", FSP, "_env_scales.csv")))

# Scale each layer in r_env by the corresponding mean/sd in scales_df
r_env_scaled <- r_env
for (v in env_vars) {
    mean_val <- scales_df$mean[scales_df$variable_name == v]
    sd_val <- scales_df$sd[scales_df$variable_name == v]
    if (v == "Intercept") {
        r_env_scaled[[v]] <- 1
    } else {
        r_env_scaled[[v]] <- (r_env[[v]] - mean_val) / sd_val
    }
}
r_env <- r_env_scaled
rm(r_env_scaled)

# ---- Effort Data ----
r_eff <- r[[c("num_observers", "duration", "distance")]]

# Scale effort layers using means/sds from scales_df
effort_vars <- c("num_observers", "duration", "distance")
r_eff_scaled <- r_eff
for (v in effort_vars) {
    mean_val <- scales_df$mean[scales_df$variable_name == v]
    sd_val <- scales_df$sd[scales_df$variable_name == v]
    r_eff_scaled[[v]] <- (r_eff[[v]] - mean_val) / sd_val
}
r_eff_log <- log(abs(r_eff_scaled[[1]] * r_eff_scaled[[2]] * r_eff_scaled[[3]]))

area_vect <- terra::vect(file.path(raw_directory, "area.shp"))

# ---- Main Loop ----
TAGS <- c("PHYLO_BASE", "PHYLO_ALL", "NOPHYLO_BASE")
for (tag in TAGS) {
    MODEL_NAME <- if (tag == "NOPHYLO_BASE") "LGCP_nophylo_background" else "LGCP_background"
    print(tag)
    print("Running 033-spatial_prediction.R ....")

    shortname_model <- paste0(data_name, "_", MODEL_NAME)
    res_dir <- file.path(res_directory, shortname)
    load(file.path(root, "res", shortname, paste0(shortname_model, ".Rdata"))) # loads 'all'

    fsps <- if (FSP == "ALL") spList[[CLUSTER]] else FSP
    failed_species <- character(0)
    # seq_along(fsps)
    for (i in 1) {
        fsp <- fsps[i]
        print(sprintf("Running for %s ...%d/%d", fsp, i, length(fsps)))
        tryCatch(
            {
                # ---- Expert Range ----
                range <- terra::vect(file.path(expert_directory, fsp, paste0(fsp, ".shp")))

                # ------ Data for extent -----
                contents <- load(file.path(
                    raw_directory, "hummingbirdSA", CLUSTER,
                    paste0(DEF_LEV, "_", FSP, "_run_files.Rdata")
                )) # loads obj store

                species_index <- which(store$species_index_map == fsp)
                sites_species <- store$y$site[store$y$species == species_index]
                cood_species <- store$cood[store$cood$site %in% sites_species, c("lon", "lat")]
                extent <- getExtentDf(cood_species, c("lon", "lat"))
                extent <- terra::ext(extent)
                # cood_vect <- terra::vect(cood_species, crs = "EPSG:4326")
                # conv_hull <- terra::convHull(cood_vect)

                # ---- Posterior Extraction ----
                stan_fit <- all$fit
                posterior_B <- rstan::extract(stan_fit, pars = "B")$B
                posterior_sigma_f <- rstan::extract(stan_fit, pars = "sigma_f")$sigma_f
                num_samples <- dim(posterior_B)[1]
                num_species <- dim(posterior_B)[2]
                num_predictors <- dim(posterior_B)[3]

                # ---- Conditional Predictions ----
                if (tag == "PHYLO_ALL") {
                    cond_filename <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO, "_cond_pred.Rdata")
                    load(file.path(cond_pred_directory, cond_filename)) # loads B_cond_list
                    num_samples <- dim(posterior_B)[1]
                    num_species <- dim(posterior_B)[2]
                    num_predictors <- dim(posterior_B)[3]
                    B_cond_array <- array(NA, dim = c(num_samples, num_species, num_predictors))
                    for (j in seq_len(num_species)) {
                        B_cond_array[, j, ] <- as.matrix(B_cond_list[[j]])
                    }
                    posterior_B <- B_cond_array
                }
                means_B <- colMeans(posterior_B[, i, ])

                if (SOFT_CLIP_IN_PRED) {
                    soft_clip <- terra::rast(file.path("/Volumes/LaCie/phylo-sdms2/analysis", "soft_clips", paste0(fsp, "_soft_clipped.tif")))
                    soft_clip_log <- log(soft_clip + 0.0001)
                }

                # # Here is my experiment
                # lin_pred <- c(r_env_crop, r_eff_scaled, soft_clip_log)

                # # ---- Linear Predictor ----

                # # posterior over B and f, but streamed (no stacks)
                # S <- min(10, dim(posterior_B)[1])
                # draws <- round(seq(1, dim(posterior_B)[1], length.out = S))
                # posterior_B_thin <- posterior_B[draws, , , drop = FALSE]

                # # exact marginal over f for each draw: exp(0.5 * sigma_f^2)
                # c_f <- exp(0.5 * posterior_sigma_f^2)

                # mean_lambda <- NULL

                # for (k in seq_along(draws)) {
                #     s <- draws[k]
                #     beta_s <- posterior_B[s, i, ] # species i
                #     # Xβ via terra (vectorized across layers)
                #     linpred_s <- terra::app(lin_pred, function(v) as.numeric(crossprod(v, beta_s)),
                #         cores = 1
                #     ) # cores>1 is fine too if you have RAM
                #     # add offset
                #     linpred_s <- linpred_s + r_eff_log_crop + soft_clip_log
                #     # exact expectation over f_s (no random draws)
                #     lambda_s <- terra::app(linpred_s, function(x) exp(x) * c_f[s])

                #     if (is.null(mean_lambda)) {
                #         mean_lambda <- lambda_s
                #     } else {
                #         mean_lambda <- ((k - 1) / k) * mean_lambda + (1 / k) * lambda_s
                #     }
                #     rm(linpred_s, lambda_s)
                #     gc()
                # }
                # print(Sys.time())
                # mean_lambda is your posterior mean intensity

                # # ex <- terra::ext(range)
                # # extent <- terra::ext(ex[1] - 4, ex[2] + 4, ex[3] - 4, ex[4] + 4)
                # # xmin <- extent[1]; xmax <- extent[2]; ymin <- extent[3]; ymax <- extent[4]
                # # r_env_crop <- terra::crop(r_env, extent)
                r_env_crop <- terra::mask(r_env, range)
                # # r_eff_log_crop <- terra::crop(r_eff_log, extent)
                r_eff_log_crop <- terra::mask(r_eff_log, range)
                r_eff_scaled <- terra::mask(r_eff_scaled, range)

                log_lambda_raster_crop <- sum(r_env_crop * means_B)

                # # f_new <- rnorm(1, mean = 0, sd = mean(posterior_sigma_f))
                if (SOFT_CLIP_IN_PRED) {
                    soft_clip_log_crop <- terra::mask(soft_clip_log, range)
                    log_lambda_raster <- log_lambda_raster_crop + soft_clip_log_crop + r_eff_log_crop

                    # Doing it the way predict_fn does it
                    #     log_lambda_raster_samples <- lapply(1:num_samples, function(s) {
                    #         sum(r_env_crop * posterior_B[s, i, ]) +
                    #             rnorm(1, 0, posterior_sigma_f[s]) + r_eff_log_crop + soft_clip_log_crop
                    #     })
                    #     log_lambda_raster <- mean(stack(log_lambda_raster_samples))
                    # } else {
                    # log_lambda_raster <- log_lambda_raster_crop + r_eff_log_crop
                    # Doing it the way predict_fn does it
                    # log_lambda_raster_samples <- lapply(1:num_samples, function(s) {
                    #     sum(r_env_crop * posterior_B[s, i, ]) +
                    #         rnorm(1, 0, posterior_sigma_f[s]) + r_eff_log_crop
                    # })
                    # log_lambda_raster <- mean(stack(log_lambda_raster_samples))
                }

                # ---- Intensity & Relative Probability ----
                # log_lambda_raster <- mean_lambda
                lambda_raster <- exp(log_lambda_raster)
                relative_prob <- lambda_raster / terra::global(lambda_raster, max, na.rm = TRUE)[1, 1]
                corrected <- relative_prob

                # ---- Output Filenames ----
                base_fn <- paste(exp_name, CLUSTER, fsp, DEF_LEV, REPNO, tag, sep = "_")
                corrected_fn <- paste0(base_fn, "_NOT_CORRECTED.tif")
                terra::writeRaster(corrected, file.path(output_dir, corrected_fn), overwrite = TRUE)
                terra::plot(corrected, main = paste("Relative Probability -", fsp, "-", tag, "-", DEF_LEV))

                # ---- Hard Clip ----
                if (HARD_CLIP) {
                    hard_clip_fn <- paste0(base_fn, "_NOT_CORRECTED_hard_clip.tif")
                    hard_clip <- terra::mask(corrected, range)
                    terra::writeRaster(hard_clip, file.path(output_dir, hard_clip_fn), overwrite = TRUE)
                    if (BINARY_CLIP) {
                        binary_fn <- paste0(base_fn, "_NOT_CORRECTED_BINARY_hard_clip.tif")
                        binary_raster <- hard_clip > 0
                        terra::writeRaster(binary_raster, file.path(output_dir, binary_fn), overwrite = TRUE)
                    }
                }

                # ---- Binary Non-Clip ----
                if (BINARY_NONCLIP) {
                    binary_fn <- paste0(base_fn, "_NOT_CORRECTED_BINARY_non_clip.tif")
                    binary_raster <- corrected > 0
                    terra::writeRaster(binary_raster, file.path(output_dir, binary_fn), overwrite = TRUE)
                }

                if (SOFT_CLIP) {
                    soft_clip <- terra::rast(file.path("/Volumes/LaCie/phylo-sdms2/analysis", "soft_clips", paste0(fsp, "_soft_clipped.tif")))
                    corrected_soft <- corrected * soft_clip
                    terra::writeRaster(corrected_soft, file.path(output_dir, paste0(base_fn, "_soft_clipped.tif")), overwrite = TRUE)
                }

                rm(corrected, relative_prob, lambda_raster, log_lambda_raster, binary_raster, hard_clip)
                gc()
            },
            error = function(e) {
                message(sprintf("Error processing species: %s : %s", fsp, e$message))
                failed_species <<- c(failed_species, fsp)
            }
        )
    }

    if (length(failed_species) > 0) {
        message("The following species failed:")
        print(failed_species)
    }
}
