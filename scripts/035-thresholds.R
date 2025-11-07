# This script computes per-species thresholds (and AUC, sens, spec)
# from existing spatial prediction rasters to support binary maps.

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(terra)
    library(pROC)
})

# ----------------------------- Configuration ------------------------------

HPC <- Sys.getenv("HPC")
DATA_SAVER <- TRUE
EXP_ROOT <- "happy"

# Roots and output locations (local vs HPC)
if (HPC == "FALSE") {
    root <- "~/phylo-sdms2"
    message("Running locally")
    tempdir <- function() "/Volumes/LaCie/temp"
    output_dir <- file.path("/Volumes/LaCie/phylo-sdms2", "analysis", EXP_ROOT, "spatial_pred")
} else {
    root <- "/vast/palmer/pi/jetz/ss4224/phylo-sdms2"
    message("Running on HPC")
    tempdir <- function() "/vast/palmer/scratch/jetz/ss4224"
    output_dir <- "/vast/palmer/scratch/jetz/ss4224/spatial_pred"
}

# Prefer external tempdir if present (local)
if (dir.exists("/Volumes/LaCie/temp")) {
    terra::terraOptions(tempdir = "/Volumes/LaCie/temp")
} else {
    message("External temp drive not found. Using default tempdir.")
}

message("tempdir(): ", tempdir())
message("output_dir: ", output_dir)

# Handy cleaner (optional)
cleanup_terra_temp <- function(path = terra::terraOptions()$tempdir, days_old = 0) {
    files <- list.files(path, full.names = TRUE)
    if (!length(files)) {
        return(invisible())
    }
    old_files <- files[file.info(files)$mtime < (Sys.time() - days_old * 86400)]
    if (length(old_files)) unlink(old_files, recursive = TRUE)
    message(length(old_files), " files removed from terra tempdir.")
}

# ----------------------------- Paths --------------------------------------

scripts_directory <- file.path(root, "scripts")
data_directory <- file.path(root, "data")
raw_directory <- file.path(root, "raw_data")
res_directory <- file.path(root, "res")
analysis_directory <- file.path(root, "analysis")
expert_directory <- file.path(root, "expert_ranges")

pred_rast_directory <- if (DATA_SAVER) {
    if (HPC == "FALSE" || HPC == "") {
        file.path("/Volumes/LaCie/phylo-sdms2", "analysis", "pred_rast")
    } else {
        file.path(root, "analysis", "pred_rast")
    }
} else {
    file.path(analysis_directory, "pred_rast")
}

# ------------------------ Experiment parameters ---------------------------

if (interactive()) {
    EXP_ROOT <- "happy" # Root name for experiment
    EXP_ID <- "woopsie"
    FSP <- "ALL"
    DEF_LEV <- "e2024"
    REPNO <- "1"
    CLUSTER <- "Coeligena" # one cluster name or NULL for all
    TEMPORAL <- FALSE
    ART <- FALSE
} else {
    args <- commandArgs(trailingOnly = TRUE)
    EXP_ROOT <- args[1]
    EXP_ID <- args[2]
    FSP <- args[3]
    DEF_LEV <- args[4]
    REPNO <- args[5]
    CLUSTER <- args[6]
    TEMPORAL <- as.logical(args[7])
    ART <- as.logical(args[8])
}
DEBUG <- FALSE
if (DEBUG) {
    print("DEBUG mode is ON")
    output_dir <- "~/Downloads/flash"
}

# Derived identifiers (per-cluster later)
exp_name <- paste0(EXP_ROOT, "_", EXP_ID)

# ------------------------------ Helpers -----------------------------------

# Decide model label by TAG (string → filename part)
model_name_for_tag <- function(tag) {
    if (tag == "NOPHYLO_BASE") "LGCP_nophylo_background" else "LGCP_background"
}

# Prediction raster path helper
pred_file_path <- function(out_dir, shortname, species, def_lev, repno, tag, trailing = "relprob") {
    file.path(out_dir, paste0(shortname, "_", species, "_", def_lev, "_", repno, "_", tag, "_", trailing, ".tif"))
}

# Build species test rows and coordinates for extraction
build_species_test <- function(store, list_of_indices, repno, species_idx) {
    idx_test <- sort(list_of_indices[[as.numeric(repno)]]$testing_sites)
    y <- store$y
    cood <- store$cood

    y_test <- y[y$site %in% idx_test, , drop = FALSE]
    rows <- which(y_test$species == species_idx)
    if (!length(rows)) {
        return(NULL)
    }

    list(
        test_binary = as.integer(y_test$count[rows] > 0),
        cood_test   = cood[y_test$site, , drop = FALSE][rows, c("lon", "lat"), drop = FALSE]
    )
}

# Compute ROC/AUC and best-threshold metrics
roc_metrics <- function(obs_binary, pred_numeric) {
    roc <- pROC::roc(response = obs_binary, predictor = pred_numeric)
    auc_val <- unlist(pROC::auc(roc))
    sens_vals <- unlist(pROC::coords(roc, "best", ret = "sensitivity"))
    spec_vals <- unlist(pROC::coords(roc, "best", ret = "specificity"))
    thr_vals <- unlist(pROC::coords(roc, "best", ret = "threshold"))
    list(
        auc  = if (length(auc_val) > 1) max(auc_val, na.rm = TRUE) else auc_val,
        sens = if (length(sens_vals) > 1) max(sens_vals, na.rm = TRUE) else sens_vals,
        spec = if (length(spec_vals) > 1) max(spec_vals, na.rm = TRUE) else spec_vals,
        thr  = if (length(thr_vals) > 1) max(thr_vals, na.rm = TRUE) else thr_vals
    )
}

# ------------------------------ Load base data -----------------------------

source(file.path(scripts_directory, "000-phyloGenie_functions.R"))
source(file.path(scripts_directory, "032-eval_functions.R"))

# species list by cluster
load(file.path(raw_directory, "hummingbirdSA", "spList.Rdata")) # loads spList

cluster_list <- if (CLUSTER == "NULL" || is.null(CLUSTER)) names(spList) else CLUSTER

TAGS <- c("PHYLO_BASE", "PHYLO_ALL", "NOPHYLO_BASE")

# ================================ MAIN =====================================

for (cluster in cluster_list) {
    tryCatch(
        {
            message("Processing cluster: ", cluster)
            fspALL <- spList[[cluster]]
            J <- length(fspALL)

            shortname <- paste0(exp_name, "_", cluster)
            data_name <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO)

            # Load store + indices for this cluster
            if (TEMPORAL) {
                stop("TEMPORAL path requires TEMPORAL_TEST_DEF_LEV; omitted here to preserve behavior.")
            } else {
                load(file.path(raw_directory, "hummingbirdSA", cluster, paste0(DEF_LEV, "_", FSP, "_run_files.Rdata"))) # store
            }
            load(file.path(raw_directory, "hummingbirdSA", cluster, paste0(DEF_LEV, "_", FSP, "_indices.Rdata"))) # list_of_indices

            species_indices <- if (ART) {
                unique(store$y$species[store$y$site %in% list_of_indices[[as.numeric(REPNO)]]$testing_sites])
            } else {
                seq_len(J)
            }

            for (tag in TAGS) {
                message(sprintf("Cluster %s | TAG %s", cluster, tag))

                # containers
                auc_values <- rep(NA_real_, J)
                sensitivity <- rep(NA_real_, J)
                specificity <- rep(NA_real_, J)
                threshold <- rep(NA_real_, J)

                for (species_idx in species_indices) {
                    fsp <- fspALL[species_idx]
                    print(fsp)

                    # Build test rows for this species
                    st <- build_species_test(store, list_of_indices, REPNO, species_idx)
                    if (is.null(st)) next

                    # Locate prediction raster and extract values at test points
                    pred_file <- pred_file_path(output_dir, shortname, fsp, DEF_LEV, REPNO, tag, trailing = "relprob")
                    if (!file.exists(pred_file)) {
                        warning("Missing spatial prediction raster: ", pred_file)
                        next
                    }
                    spr <- terra::rast(pred_file)
                    vals <- terra::extract(spr, st$cood_test, method = "simple", ID = FALSE)

                    # Safety: prediction column name
                    pred_num <- suppressWarnings(as.numeric(vals[[1]]))
                    if (!length(pred_num)) {
                        warning("Predictions missing/NA for ", fsp, " at ", tag)
                        next
                    }

                    # Metrics
                    m <- roc_metrics(st$test_binary, pred_num)
                    auc_values[species_idx] <- m$auc
                    sensitivity[species_idx] <- m$sens
                    specificity[species_idx] <- m$spec
                    threshold[species_idx] <- m$thr

                    # Other thresholds
                    # target_prev <- mean(st$test_binary == 1) # or your independent prevalence estimate
                    # thr_prev <- as.numeric(quantile(pred_num, probs = 1 - target_prev, na.rm = TRUE))
                    # thr_p10 <- as.numeric(quantile(pred_num[st$test_binary == 1], probs = 0.050, na.rm = TRUE))
                    # threshold[species_idx] <- thr_p10

                    # Write thresholds table for this cluster + tag
                    threshold_dir <- file.path(analysis_directory, EXP_ROOT, "eval", "thresholds")
                    dir.create(threshold_dir, recursive = TRUE, showWarnings = FALSE)

                    out_csv <- file.path(
                        threshold_dir,
                        paste0(shortname, "_", fsp, "_", DEF_LEV, "_", REPNO, "_", tag, "_thresholds.csv")
                    )

                    write.csv(data.frame(
                        species     = fspALL,
                        AUC         = auc_values,
                        sensitivity = sensitivity,
                        specificity = specificity,
                        threshold   = threshold
                    ), out_csv, row.names = FALSE)

                    message("Wrote: ", out_csv)
                }
            }
        },
        error = function(e) {
            message("Error occurred in cluster ", cluster, ": ", conditionMessage(e))
            traceback()
        }
    )
}
